#ifndef __STEREO_COMPUTE_H__
#define __STEREO_COMPUTE_H__

#include <fstream>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"

#include "mio/altro/error.h"
#include "mio/altro/casting.h"
#include "mio/altro/macros.h"
#include "mio/altro/timers.h"
#include "mio/altro/types.h"
#include "mio/math/math.h"
#include "mvg/stereo_types.h"
#include "mvg/stereo_io.h"


template<typename TYPE>
inline void ComputeEssentialFundamentalT(const cv::Mat &R, const cv::Mat &T, const cv::Mat &M0, const cv::Mat &M1, 
                                         cv::Mat &E, cv::Mat &F){
  TYPE *T_ = (TYPE *)(T.data);
  TYPE S_[] = {0, -T_[2], T_[1], T_[2], 0, -T_[0], -T_[1], T_[0], 0};
  cv::Mat S(3, 3, R.type(), S_);
  E = S*R;
  F = M1.inv().t() * E * M0.inv();
  TYPE val = F.at<TYPE>(8);
  if(val > 0)
    F *= 1.0/val;
}

//given R, T, cam_matrix[0], and cam_matrix[1], this function computes the essential (E) and fundamental (F) matrices
inline void ComputeEssentialFundamental(stereoCalData_t &cal_data){
  const int depth = cal_data.R.depth();
  if(depth == CV_32F)
    ComputeEssentialFundamentalT<float>(cal_data.R, cal_data.T, cal_data.K[0], 
                                        cal_data.K[1], cal_data.E, cal_data.F);
  else if(depth == CV_64F)
    ComputeEssentialFundamentalT<double>(cal_data.R, cal_data.T, cal_data.K[0], 
                                         cal_data.K[1], cal_data.E, cal_data.F);
}


//which_image = 0: pnt is in left image, computed line exists in right image
//which_image = 1: pnt is in right image, computed line exists in left image
template <typename PNT_T, typename LINE_T>
void FindEpipolarLine(const PNT_T &pnt, const cv::Mat &F, LINE_T &line, const int which_image){
  STD_INVALID_ARG_E(which_image == 0 || which_image == 1)
  const int F_depth = F.depth();
  STD_INVALID_ARG_E(F_depth == CV_32F || F_depth == CV_64F)
#define FIND_EPIPOLAR_LINE(data_type)\
    data_type *fd = mio::StaticCastPtr<data_type>(F.data);\
    if(which_image == 0){ /* F * pnt */ \
      line.a = fd[0]*pnt.x + fd[1]*pnt.y + fd[2];\
      line.b = fd[3]*pnt.x + fd[4]*pnt.y + fd[5];\
      line.c = fd[6]*pnt.x + fd[7]*pnt.y + fd[8];\
    }\
    else{ /* F.t() * pnt */ \
      line.a = fd[0]*pnt.x + fd[3]*pnt.y + fd[6];\
      line.b = fd[1]*pnt.x + fd[4]*pnt.y + fd[7];\
      line.c = fd[2]*pnt.x + fd[5]*pnt.y + fd[8];\
    }
  if(F_depth == CV_32F){
    FIND_EPIPOLAR_LINE(float)
  }
  else{
    FIND_EPIPOLAR_LINE(double)
  }
#undef FIND_EPIPOLAR_LINE
}

template <typename PNT_T>
inline void DrawEpipolarLine(cv::Mat &img, const cv::Mat &F, const PNT_T &pnt, const int which_image){
  STD_INVALID_ARG_E(which_image == 0 || which_image == 1)
  line3d_t line;
  FindEpipolarLine(pnt, F, line, which_image);
  cv::line( img, sm::SolveLineEqn2<cv::Point>(line, 0),
            sm::SolveLineEqn2<cv::Point>(line, img.cols), cv::Scalar::all(255) );
}


inline void ComputeProjectionMatrices(const stereoCalData_t &cal_data, cv::Mat &P1, cv::Mat &P2){
  STD_INVALID_ARG_E(cal_data.K[0].rows == 3 && cal_data.K[0].cols == 3)
  STD_INVALID_ARG_E(cal_data.K[1].rows == 3 && cal_data.K[1].cols == 3)
  STD_INVALID_ARG_E(cal_data.T.rows == 3 && cal_data.T.cols == 1)

  cv::Mat M_ext = cv::Mat::eye(3, 4, CV_64F);
  P1 = cal_data.K[0] * M_ext;
  cv::Mat R_roi( M_ext, cv::Rect( cv::Point(0, 0), cv::Size(3, 3) ) ),
          T_roi( M_ext, cv::Rect( cv::Point(3, 0), cv::Size(1, 3) ) );
  cal_data.R.copyTo(R_roi);
  cal_data.T.copyTo(T_roi);
  P2 = cal_data.K[1] * M_ext;
}


inline void TriangulatePoints(const cv::Mat &P1, const cv::Mat &P2,
                              const std::vector<cv::Point2f> &matched_pixels1,
                              const std::vector<cv::Point2f> &matched_pixels2, cv::Mat &points4D){
  STD_INVALID_ARG_E(P1.depth() == CV_64F && P2.depth() == CV_64F)
  STD_INVALID_ARG_E(P1.rows == 3 && P1.cols == 4 && P2.rows == 3 && P2.cols == 4)
  STD_INVALID_ARG_E( matched_pixels1.size() > 0 && matched_pixels1.size() == matched_pixels2.size() )

  const size_t num_match = matched_pixels1.size();
  const double *P1_data = mio::StaticCastPtr<double>(P1.data),
               *P2_data = mio::StaticCastPtr<double>(P2.data);
  points4D = cv::Mat(num_match, 4, CV_32F);

  auto tbb_func = [&](const tbb::blocked_range<size_t> &br){
    //'A' is the cross product of the image point [x; y; 1] and (PX) for P1 and P2
    cv::Mat A(6, 4, CV_64F), //m×n real or complex matrix
            D, //m×n rectangular diagonal matrix with non-negative real numbers. diagonal is the singular values of A.
            U, //m×m real or complex unitary matrix. columns of U are left-singular vectors of A.
            Vt; //n×n real or complex unitary matrix. columns of V are the right-singular vectors of A.
    double *a = mio::StaticCastPtr<double>(A.data);

    const size_t br_end = br.end();
    for(size_t i = br.begin(); i != br_end; ++i){
      {
        const double x = matched_pixels1[i].x,
                     y = matched_pixels1[i].y;
        const double *p = P1_data;
        a[0] = y*p[8] -   p[4];  a[1] = y*p[9] -   p[5];  a[2] =  y*p[10] -  p[6];   a[3] =  y*p[11] -  p[7];
        a[4] =   p[0] - x*p[8];  a[5] =   p[1] - x*p[9];  a[6] =    p[2] - x*p[10];  a[7] =    p[3] - x*p[11];
        a[8] = x*p[4] - y*p[0];  a[9] = x*p[5] - y*p[1];  a[10] = x*p[6] - y*p[2];   a[11] = x*p[7] - y*p[3];
      }
      const double x = matched_pixels2[i].x,
                   y = matched_pixels2[i].y;
      const double *p = P2_data;
      a[12] = y*p[8] -   p[4];  a[13] = y*p[9] -   p[5];  a[14] = y*p[10] -  p[6];   a[15] = y*p[11] -  p[7];
      a[16] =   p[0] - x*p[8];  a[17] =   p[1] - x*p[9];  a[18] =   p[2] - x*p[10];  a[19] =   p[3] - x*p[11];
      a[20] = x*p[4] - y*p[0];  a[21] = x*p[5] - y*p[1];  a[22] = x*p[6] - y*p[2];   a[23] = x*p[7] - y*p[3];

      //the null space of A is spanned by the last n − r columns of V (n is square dim. of V, r is rank of V)
      //we have Vt so take last row to get null space of A
      cv::SVD::compute(A, D, U, Vt, cv::SVD::MODIFY_A);
      Vt.row(3).copyTo( points4D.row(i) );
    }
  };

  tbb::parallel_for(tbb::blocked_range<size_t>(0, num_match), tbb_func);
}


inline void TriangulatePoints(const cv::Mat &P1, const cv::Mat &P2,
                              const std::vector<cv::Point2f> &matched_pixels1,
                              const std::vector<cv::Point2f> &matched_pixels2,
                              std::vector<vertex3f_t> &out_pnts){
  STD_INVALID_ARG_E(P1.rows == 3 && P1.cols == 4 && P2.rows == 3 && P2.cols == 4)
  STD_INVALID_ARG_E( matched_pixels1.size() > 0 && matched_pixels1.size() == matched_pixels2.size() )
  cv::Mat points4D;
  //cv::triangulatePoints(P1, P2, matched_pixels1, matched_pixels2, points4D); //points4D is 4xN
  //OpenCV docs says .t() does not perform an actual transpose, the assignment operator will take care of that
  //const cv::Mat points4D_t = points4D.t();
  TriangulatePoints(P1, P2, matched_pixels1, matched_pixels2, points4D); //points4D is Nx4
  STD_INVALID_ARG_E(points4D.depth() == CV_32F)

  const size_t num_pnts = points4D.rows;
  const vertex4f_t *pnt4D = mio::StaticCastPtr<vertex4f_t>(points4D.data);

  out_pnts.resize(num_pnts);
  for(size_t i = 0; i < num_pnts; ++i){
    const float scale = (pnt4D[i].w != 0.0f) ? 1.0f / pnt4D[i].w : 1.0f;
    out_pnts[i] = vertex3f_t(pnt4D[i].x*scale, pnt4D[i].y*scale, pnt4D[i].z*scale);
  }
}


inline int ExtractCalTargetPoints(const std::string &file_path, //in
                                  const camCalTarget_t &cal_target, //in
                                  const std::vector<std::string> &img_file_name_vec, //in
                                  stereoCalData_t &cal_data, //out
                                  const int find_target_flags,
                                  const size_t num_camera){
  STD_INVALID_ARG_E(num_camera == 1 || num_camera == 2)
  STD_INVALID_ARG_E(cal_target.type_str == "chess" || cal_target.type_str == "circle" || cal_target.type_str == "a-circle")
  const int max_scale = 3;
  const bool display_corners = false;
  const size_t num_img_per_cam = img_file_name_vec.size()/num_camera;
  size_t num_valid_img_pair = 0;

  cal_data.img_cal_pnts[0].resize(num_img_per_cam);
  cal_data.img_cal_pnts[1].resize(num_img_per_cam);
  cal_data.good_img_file_names[0].resize(num_img_per_cam);
  cal_data.good_img_file_names[1].resize(num_img_per_cam);
  cal_data.img_size = cv::Size();

  //extract calibration target features
  for(size_t i = 0, j; i < num_img_per_cam; ++i){
    for(j = 0; j < num_camera; ++j){
      const std::string &file_name = img_file_name_vec[i*num_camera + j];
      cv::Mat cal_img = LoadCalImgDefaultProc(file_path + "/" + file_name);
      if( cal_img.empty() )
        break;
      if( cal_data.img_size == cv::Size() )
        cal_data.img_size = cal_img.size(); //initialize cal_data.img_size
      else if(cal_img.size() != cal_data.img_size){
        std::cout << "The image " << file_name << " has a different size than the first image. skipping image pair.\n";
        break;
      }

      bool found = false;
      std::vector<cv::Point2f> &img_cal_pnts = cal_data.img_cal_pnts[j][num_valid_img_pair]; //make a reference to the corner vector for this image
      const size_t scale_start = (cal_img.cols < 400) ?  2 : 1;
      for(size_t scale = scale_start; scale <= max_scale; scale++){
        cv::Mat cal_img_scaled;
        if(scale == 1)
          cal_img_scaled = cal_img;
        else
          cv::resize(cal_img, cal_img_scaled, cv::Size(), scale, scale);

        if(cal_target.type_str == "chess")
          found = cv::findChessboardCorners(cal_img_scaled, cal_target.size, img_cal_pnts, find_target_flags);
        else
          found = cv::findCirclesGrid(cal_img_scaled, cal_target.size, img_cal_pnts, find_target_flags);

        if(found){
          if(scale > 1){
            cv::Mat cornersMat(img_cal_pnts);
            cornersMat *= 1./scale;
          }
          break;
        }
      }
      if(!found){
        std::cout << "No calibration target found in " << file_name << ", skipping image pair" << std::endl;
        break;
      }

      if(cal_target.type_str == "chess")
        cv::cornerSubPix( cal_img, img_cal_pnts, cv::Size(11, 11), cv::Size(-1, -1), 
                          cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 30, 0.01) );
    }
    if(j == num_camera){
      for(size_t k = 0; k < num_camera; ++k)
        cal_data.good_img_file_names[k][num_valid_img_pair] = img_file_name_vec[i*num_camera + k];
      ++num_valid_img_pair;
    }
  }

  EXP_CHK( num_valid_img_pair > 1, return(-1) )

  cal_data.img_cal_pnts[0].resize(num_valid_img_pair);
  cal_data.img_cal_pnts[1].resize(num_valid_img_pair);
  cal_data.good_img_file_names[0].resize(num_valid_img_pair);
  cal_data.good_img_file_names[1].resize(num_valid_img_pair);
  cal_data.SetNumValidImgPair(num_valid_img_pair);

  return num_valid_img_pair;
}


inline void CalibrateCamera(const camCalTarget_t &cal_target, stereoCalData_t &cal_data,
                            const size_t cam_idx, int cal_flags, int *rms_ = nullptr,
                            std::vector<float> *reproj_errs_ = nullptr){
  STD_INVALID_ARG_E(cam_idx == 0 || cam_idx == 1)
  std::vector<cv::Mat> rvecs, tvecs;
  CHK_PTR_MAKE_REF(int, rms_, rms)
  CHK_PTR_MAKE_REF(std::vector<float>, reproj_errs_, reproj_errs)

  cal_flags &= ~cv::CALIB_USE_INTRINSIC_GUESS;
  //these options are only supported by cv::stereoCompute()
  cal_flags &= ~cv::CALIB_FIX_INTRINSIC;
  cal_flags &= ~cv::CALIB_FIX_FOCAL_LENGTH;
  cal_flags &= ~cv::CALIB_SAME_FOCAL_LENGTH;

  std::vector< std::vector<cv::Point3f> > target_points( cal_data.GetNumImgPerCam() );
  std::fill(target_points.begin(), target_points.end(), cal_target.point_vec);

  rms = cv::calibrateCamera(target_points, cal_data.img_cal_pnts[cam_idx], cal_data.img_size,
                            cal_data.K[cam_idx], cal_data.D[cam_idx], rvecs, tvecs,
#if ICV_OPENCV_VERSION_MAJOR < 3
                            cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, 100, 1e-5),
                            cal_flags);
#else 
                            cal_flags,
                            cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, 100, 1e-5) );
#endif

//  double totalAvgErr = cv::computeReprojectionErrors(target_points, cal_data.img_cal_pnts[cam_idx], rvecs, tvecs,
//                                                     cal_data.K[cam_idx], cal_data.D[cam_idx], reproj_errs);
  CHK_PTR_MAKE_REF_CLEAN_UP(rms_)
  CHK_PTR_MAKE_REF_CLEAN_UP(reproj_errs_)
}


inline void StereoCalibrate(const camCalTarget_t &cal_target,
                            stereoCalData_t &cal_data,
                            std::vector<double> &repro_err_vec, double &repro_err_1, double &repro_err_2,
                            const int cal_flags, const bool pre_cal_cameras){
  std::cout << "Running stereo calibration ...\n";
    
  std::vector< std::vector<cv::Point3f> > target_points( cal_data.GetNumImgPerCam() );
  std::fill(target_points.begin(), target_points.end(), cal_target.point_vec);

  if(pre_cal_cameras){
    CalibrateCamera(cal_target, cal_data, 0, cal_flags);
    CalibrateCamera(cal_target, cal_data, 1, cal_flags);
    mio::Print(cal_data.K[0], "Pre-M1");
    mio::Print(cal_data.K[1], "Pre-M2");
  }

  repro_err_1 = cv::stereoCalibrate( target_points, cal_data.img_cal_pnts[0], cal_data.img_cal_pnts[1],
                                     cal_data.K[0], //out - intrinsic parameters camera1
                                     cal_data.D[0],   //out - distortion coefficients camera 1
                                     cal_data.K[1], //out - intrinsic parameters camera2
                                     cal_data.D[1],   //out - distortion coefficients camera 1
                                     cal_data.img_size, 
                                     cal_data.R, //out - rotation matrix between the 1st and the 2nd camera coordinate systems
                                     cal_data.T, //out - translation vector between the coordinate systems of the camera
                                     cal_data.E, //out - essential matrix
                                     cal_data.F, //out - fundamental matrix
#if ICV_OPENCV_VERSION_MAJOR < 3
                                     cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, 100, 1e-5),
                                     cal_flags
#else 
                                     cal_flags,
                                     cv::TermCriteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, 100, 1e-5)
#endif
                                     );



  //calibration quality check
  //because the output fundamental matrix implicitly includes all the output information,
  //we can check the quality of calibration using the epipolar geometry constraint: m2^t * F * m1 = 0
  //FIXME - 20 cal pairs: 15 for cal, 5 for below.
  double error_sum = 0;
  size_t total_num_pnt = 0;
  std::vector<cv::Vec3f> lines[2];
  repro_err_vec.clear();
  repro_err_vec.resize( cal_data.GetNumImgPerCam() );
  for(size_t i = 0; i < cal_data.GetNumImgPerCam(); i++){
    const size_t num_pnt = cal_data.img_cal_pnts[0][i].size();
    cv::Mat img_pnt_mat[2];
    for(size_t k = 0; k < 2; k++){
      img_pnt_mat[k] = cv::Mat(cal_data.img_cal_pnts[k][i]);
      cv::undistortPoints(img_pnt_mat[k], img_pnt_mat[k], cal_data.K[k],
                          cal_data.D[k], cv::Mat(), cal_data.K[k]);
      cv::computeCorrespondEpilines(img_pnt_mat[k], k+1, cal_data.F, lines[k]);
    }
    double pair_error = 0;
    for(size_t j = 0; j < num_pnt; j++)
      pair_error += fabs(cal_data.img_cal_pnts[0][i][j].x * lines[1][j][0] + cal_data.img_cal_pnts[0][i][j].y * lines[1][j][1] + lines[1][j][2]) +
                    fabs(cal_data.img_cal_pnts[1][i][j].x * lines[0][j][0] + cal_data.img_cal_pnts[1][i][j].y * lines[0][j][1] + lines[0][j][2]);
    repro_err_vec[i] = pair_error/static_cast<double>(num_pnt);
    error_sum += pair_error;
    total_num_pnt += num_pnt;
  }
  repro_err_2 = error_sum/static_cast<double>(total_num_pnt);
}


inline void ComputeRectification(const std::string file_path,
                                 stereoCalData_t &cal_data,
                                 const bool use_opencv_rectification = true,
                                 int stereo_rectify_flags = cv::CALIB_ZERO_DISPARITY,
                                 const int stereo_rectify_alpha = -1){
  if(use_opencv_rectification){
    cv::stereoRectify(cal_data.K[0], cal_data.D[0],
                      cal_data.K[1], cal_data.D[1],
                      cal_data.img_size, cal_data.R, cal_data.T,
                      cal_data.R1, //out - 3x3 rectification transform (rotation matrix) for the first camera
                      cal_data.R2, //out - 3x3 rectification transform (rotation matrix) for the second camera
                      cal_data.P1, //out - 3x4 projection matrix in the new (rectified) coordinate systems for the first camera
                      cal_data.P2, //out - 3x4 projection matrix in the new (rectified) coordinate systems for the second camera
                      cal_data.Q,  //out - 4x4 disparity-to-depth mapping matrix (see reprojectImageTo3D() ).
                      stereo_rectify_flags, stereo_rectify_alpha,
                      cal_data.img_size, &cal_data.valid_roi[0], &cal_data.valid_roi[1]);
  }
  else{
    //Hartley's method: use intrinsic parameters of each camera, but compute
    //the rectification transformation directly from the fundamental matrix
    std::cout << "Hartley rectification method\n";
    std::vector<cv::Point2f> all_img_pnt_vec[2];
    for(size_t k = 0; k < 2; ++k)
      for(size_t i = 0; i < cal_data.GetNumImgPerCam(); ++i)
        std::copy( cal_data.img_cal_pnts[k][i].begin(), cal_data.img_cal_pnts[k][i].end(), std::back_inserter(all_img_pnt_vec[k]) );

    cal_data.F = cv::findFundamentalMat(cv::Mat(all_img_pnt_vec[0]), cv::Mat(all_img_pnt_vec[1]), cv::FM_RANSAC, 0, 0);
    
    cv::stereoRectifyUncalibrated(cv::Mat(all_img_pnt_vec[0]), cv::Mat(all_img_pnt_vec[1]), 
                                  cal_data.F, cal_data.img_size, cal_data.H1, cal_data.H2, 3);

    cal_data.R1 = cal_data.K[0].inv() * cal_data.H1 * cal_data.K[0];
    cal_data.R2 = cal_data.K[1].inv() * cal_data.H2 * cal_data.K[1];
    cal_data.P1 = cal_data.K[0];
    cal_data.P2 = cal_data.K[1];
  }
}


inline void RectifyImages(const std::string file_path,
                          const stereoCalData_t &cal_data,
                          std::vector<cv::Mat> &rect_cal_img_vec){
  const bool use_opencv_rectification = true;
  //OpenCV can handle left-right or up-down camera arrangements
  const bool isVerticalStereo = fabs( cal_data.P2.at<double>(1, 3) ) > fabs( cal_data.P2.at<double>(0, 3) );
  printf("Rectification detected a %s stereo orientation\n", isVerticalStereo ? "vertical" : "horizontal");

  cv::Mat canvas;
  double sf; //scale factor
  size_t w, h;
  if(!isVerticalStereo){
    sf = 600./MAX(cal_data.img_size.width, cal_data.img_size.height);
    w = cvRound(cal_data.img_size.width*sf);
    h = cvRound(cal_data.img_size.height*sf);
    canvas.create(h, w*2, CV_8UC3);
  }
  else{
    sf = 300./MAX(cal_data.img_size.width, cal_data.img_size.height);
    w = cvRound(cal_data.img_size.width*sf);
    h = cvRound(cal_data.img_size.height*sf);
    canvas.create(h*2, w, CV_8UC3);
  }

  cv::Mat rectificationMap[2][2];
  //compute the undistortion and rectification transformation maps
  cv::initUndistortRectifyMap(cal_data.K[0], cal_data.D[0], cal_data.R1, cal_data.P1,
                              cal_data.img_size, CV_16SC2, rectificationMap[0][0], rectificationMap[0][1]);
  cv::initUndistortRectifyMap(cal_data.K[1], cal_data.D[1], cal_data.R2, cal_data.P2,
                              cal_data.img_size, CV_16SC2, rectificationMap[1][0], rectificationMap[1][1]);

  STD_INVALID_ARG_E( cal_data.good_img_file_names[0].size() == cal_data.good_img_file_names[0].size() )
  rect_cal_img_vec.resize( cal_data.GetNumImgPerCam() );
  for(size_t i = 0; i < cal_data.GetNumImgPerCam(); ++i){
    for(size_t k = 0; k < 2; ++k){
      cv::Mat img = mio::ReadImage(file_path + "/" + cal_data.good_img_file_names[k][i], 0), rectified_img, color_img;
      cv::remap(img, rectified_img, rectificationMap[k][0], rectificationMap[k][1], cv::INTER_LINEAR);
      cv::cvtColor(rectified_img, color_img, cv::COLOR_GRAY2BGR);
      cv::Mat canvasPart = !isVerticalStereo ? canvas( cv::Rect(w*k, 0, w, h) ) : canvas( cv::Rect(0, h*k, w, h) );
      cv::resize(color_img, canvasPart, canvasPart.size(), 0, 0, cv::INTER_AREA);
      //draw rectangle for valid pixel region (computed from stereoRectify())
      if(use_opencv_rectification){
        cv::Rect vroi(cvRound(cal_data.valid_roi[k].x*sf), cvRound(cal_data.valid_roi[k].y*sf),
                      cvRound(cal_data.valid_roi[k].width*sf), cvRound(cal_data.valid_roi[k].height*sf));
        cv::rectangle(canvasPart, vroi, cv::Scalar(0, 0, 255), 3, 8);
      }
    }

    //draw lines representing a rectified image
    if(!isVerticalStereo)
      for(int j = 0; j < canvas.rows; j += 16)
        cv::line(canvas, cv::Point(0, j), cv::Point(canvas.cols, j), cv::Scalar(0, 255, 0), 1, 8);
    else
      for(int j = 0; j < canvas.cols; j += 16)
        cv::line(canvas, cv::Point(j, 0), cv::Point(j, canvas.rows), cv::Scalar(0, 255, 0), 1, 8);

    rect_cal_img_vec[i] = canvas.clone();
  }
}


inline void DrawMatches(const cv::Mat &img1, const std::vector<cv::Point2f> &matched_pixels1, 
                        const cv::Mat &img2, const std::vector<cv::Point2f> &matched_pixels2,
                        cv::Mat &match_img){
  STD_INVALID_ARG_E( !img1.empty() && !img2.empty() )
  STD_INVALID_ARG_E( matched_pixels1.size() > 0 && matched_pixels1.size() == matched_pixels2.size() )

  const size_t num_match = matched_pixels1.size();
  std::vector<cv::KeyPoint> keypoints1(num_match), keypoints2(num_match);
  std::vector<cv::DMatch> descriptor_matches(num_match);

  for(size_t i = 0; i < num_match; ++i){
    descriptor_matches[i].queryIdx = descriptor_matches[i].trainIdx = i;
    keypoints1[i].pt = matched_pixels1[i];
    keypoints2[i].pt = matched_pixels2[i];
  }

  cv::drawMatches(img1, keypoints1, img2, keypoints2, descriptor_matches, match_img);
}


inline void MatchEpipolar(const cv::Mat &img1, const cv::Mat &img2, const cv::Mat &F,
                          std::vector<cv::Point2f> &matched_pixels1,
                          std::vector<cv::Point2f> &matched_pixels2,
                          const std::vector<cv::KeyPoint> &keypoints1,
                          const std::vector<cv::KeyPoint> &keypoints2,
                          cv::Mat &descriptors1,
                          cv::Mat &descriptors2,
                          std::vector<cv::DMatch> &descriptor_matches,
                          cv::Mat *draw_img = nullptr,
                          const std::string *matched_pix_file_full = nullptr,
                          const bool use_epipolar_constraint = true,
                          const double epipolar_constraint_dist_thresh = 0.5){
  //match the features
  cv::BFMatcher matcher(cv::NORM_L2, true);
  matcher.match(descriptors1, descriptors2, descriptor_matches);

  struct STwoPnt{
    cv::Point2f p1, p2;
    STwoPnt(const cv::Point2f &p1_, const cv::Point2f &p2_) : p1(p1_), p2(p2_) {}
    void Print(){ std::cout << p1 << ", " << p2 << "\n"; }
    inline bool Equals(const STwoPnt &p){ return (p1 == p.p1 || p2 == p.p2); }
  };
  const size_t num_match = descriptor_matches.size();
  std::vector<STwoPnt> matched_pix;
  matched_pix.reserve(num_match);

  //filter matches using the epipolar constraint
  std::vector<char> matches_mask(num_match, 1);
  if(use_epipolar_constraint){
    for(size_t i = 0; i < num_match; ++i){
      const cv::Point2d &_p1 = keypoints1[descriptor_matches[i].queryIdx].pt,
                        &_p2 = keypoints2[descriptor_matches[i].trainIdx].pt;
      line3d_t line1, line2;
      FindEpipolarLine(_p1, F, line2, 0);
      FindEpipolarLine(_p2, F, line1, 1);
      if(sm::DistToLine2(line1, _p1) > epipolar_constraint_dist_thresh ||
         sm::DistToLine2(line2, _p2) > epipolar_constraint_dist_thresh)
        matches_mask[i] = 0;
      else
        matched_pix.push_back( STwoPnt(_p1, _p2) );
    }
  }
  else
    for(size_t i = 0; i < num_match; ++i)
      matched_pix.push_back( STwoPnt(keypoints1[descriptor_matches[i].queryIdx].pt,
                                     keypoints2[descriptor_matches[i].trainIdx].pt) );

  //remove duplicate points
  const size_t num_epi_match = matched_pix.size();
  matched_pixels1.clear();
  matched_pixels1.reserve(num_epi_match);
  matched_pixels2.clear();
  matched_pixels2.reserve(num_epi_match);
  std::vector<bool> epi_matches_mask(num_epi_match, true);
  for(size_t i = 0; i < num_epi_match; ++i)
    if( epi_matches_mask[i] ){
      for(size_t j = i+1; j < num_epi_match; ++j)
        if( epi_matches_mask[j] && matched_pix[i].Equals( matched_pix[j] ) )
          epi_matches_mask[j] = false;
      matched_pixels1.push_back(matched_pix[i].p1);
      matched_pixels2.push_back(matched_pix[i].p2);
    }

//  cv::Mat img1_ = img1.clone();
//  cv::Mat img2_ = img2.clone();
//  for(auto &pnt : matched_pixels1)
//    cv::circle(img1_, pnt, 3, cv::Scalar::all(255), -1);
//  for(auto &pnt : matched_pixels2)
//    cv::circle(img2_, pnt, 3, cv::Scalar::all(255), -1);

  //draw the matches
  if(draw_img != nullptr)
    cv::drawMatches(img1, keypoints1, img2, keypoints2, descriptor_matches, *draw_img, 
                    cv::Scalar::all(-1), cv::Scalar::all(-1), matches_mask);

  if(matched_pix_file_full != nullptr){
    std::ofstream file(*matched_pix_file_full);
    if( file.is_open() ){
      for(int i = 0; i < matched_pixels1.size(); ++i)
        file << matched_pixels1[i].x << "," << matched_pixels1[i].y << "," <<
                matched_pixels2[i].x << "," << matched_pixels2[i].y << "\n";
      file.close();
    }
    else
      std::cout << "Unable to open file: " << *matched_pix_file_full << "\n";
  }
}


#define FEAT_MATCH_FUNC_NAME(func_name)\
inline void func_name(const cv::Mat &img1, const cv::Mat &img2, const cv::Mat &F,\
                      std::vector<cv::Point2f> &matched_pixels1,\
                      std::vector<cv::Point2f> &matched_pixels2,\
                      std::vector<cv::KeyPoint> *keypoints1_ = nullptr,\
                      std::vector<cv::KeyPoint> *keypoints2_ = nullptr,\
                      std::vector<cv::DMatch> *descriptor_matches_ = nullptr,\
                      cv::Mat *draw_img = nullptr,\
                      const std::string *matched_pix_file_full = nullptr,\
                      const bool use_epipolar_constraint = true,\
                      const double epipolar_constraint_dist_thresh = 0.5,

#define FEAT_MATCH_FUNC_TOP\
  CHK_PTR_MAKE_REF(std::vector<cv::KeyPoint>, keypoints1_, keypoints1)\
  CHK_PTR_MAKE_REF(std::vector<cv::KeyPoint>, keypoints2_, keypoints2)\
  CHK_PTR_MAKE_REF(std::vector<cv::DMatch>, descriptor_matches_, descriptor_matches)\
  cv::Mat descriptors1, /*queryDescriptors*/\
          descriptors2; /*trainDescriptors*/\
  const bool use_provided_keypoints = false;

#define FEAT_MATCH_FUNC_BOTTOM(feat_var)\
  feat_var(img1, cv::noArray(), keypoints1, descriptors1, use_provided_keypoints);\
  feat_var(img2, cv::noArray(), keypoints2, descriptors2, use_provided_keypoints);\
  MatchEpipolar(img1, img2, F, matched_pixels1, matched_pixels2, keypoints1, keypoints2,\
                descriptors1, descriptors2, descriptor_matches, draw_img, matched_pix_file_full,\
                use_epipolar_constraint, epipolar_constraint_dist_thresh);\
  CHK_PTR_MAKE_REF_CLEAN_UP(keypoints1_)\
  CHK_PTR_MAKE_REF_CLEAN_UP(keypoints2_)\
  CHK_PTR_MAKE_REF_CLEAN_UP(descriptor_matches_)

FEAT_MATCH_FUNC_NAME(SiftMatch) const int num_sift_features = 100.0,
                                const int sift_num_octave_layer = 3,
                                const double sift_contrast_thresh = 0.03,
                                const double sift_edge_thresh = 10.0,
                                const double sift_sigma = 1.2){
  FEAT_MATCH_FUNC_TOP
  cv::Ptr<cv::xfeatures2d::SIFT> sift;
  sift = cv::xfeatures2d::SIFT::create(num_sift_features, sift_num_octave_layer,
                                       sift_contrast_thresh, sift_edge_thresh, sift_sigma);
  FEAT_MATCH_FUNC_BOTTOM(sift->detectAndCompute)
}


FEAT_MATCH_FUNC_NAME(OrbMatch) const int num_feat = 500,
                               const float scale_factor = 1.2f,
                               const int num_levels = 8,
                               const int edge_thresh = 31,
                               const int first_level = 0, const int WTA_K = 2,
                               const int score_type = cv::ORB::HARRIS_SCORE,
                               const int patch_size = 31){
  FEAT_MATCH_FUNC_TOP
  cv::Ptr<cv::ORB> orb = cv::ORB::create(num_feat, scale_factor, num_levels, edge_thresh,
                                         first_level, WTA_K, score_type, patch_size);
  FEAT_MATCH_FUNC_BOTTOM(orb->detectAndCompute)
}


FEAT_MATCH_FUNC_NAME(BriskMatch) const int thresh = 30,
                                 const int octaves = 3,
                                 const float pattern_scale = 1.0f){
  FEAT_MATCH_FUNC_TOP
  cv::Ptr<cv::BRISK> brisk = cv::BRISK::create(thresh, octaves, pattern_scale);
  FEAT_MATCH_FUNC_BOTTOM(brisk->detectAndCompute)
}

#undef FEAT_MATCH_FUNC_NAME
#undef FEAT_MATCH_FUNC_TOP
#undef FEAT_MATCH_FUNC_BOTTOM


inline void ComputeDisparity(const cv::Mat &imgl, const cv::Mat &imgr, cv::Mat &img_disp){
  cv::Ptr<cv::StereoSGBM> sgbm;
  const int min_disparity = 0,
            num_disparity = ( (imgl.size().width/8) + 15 ) & -16,
            block_size = 3,
            num_channel = imgl.channels();
  sgbm = cv::StereoSGBM::create(min_disparity, num_disparity, block_size);
  sgbm->setP1(8 *  num_channel*block_size*block_size);
  sgbm->setP2(32 * num_channel*block_size*block_size);
  sgbm->setDisp12MaxDiff(-1);
  sgbm->setPreFilterCap(63);
  sgbm->setUniquenessRatio(10);
  sgbm->setSpeckleWindowSize(100);
  sgbm->setSpeckleRange(32);

  sgbm->setMode(cv::StereoSGBM::MODE_HH);
  sgbm->compute(imgl, imgr, img_disp);
}


#ifdef HAVE_CUDA
inline void ComputeDisparityCuda(const cv::Mat &imgl, const cv::Mat &imgr, cv::Mat &img_disp){
  cv::cuda::GpuMat d_left, d_right, d_disp(imgl.size(), CV_8U);
  d_left.upload(imgl);
  d_right.upload(imgr);
  int ndisp, iters, levels, nr_plane;
  cv::cuda::StereoConstantSpaceBP::estimateRecommendedParams(m_cal_data.img_size.width, m_cal_data.img_size.height,
                                                             ndisp, iters, levels, nr_plane);
  cv::Ptr<cv::cuda::StereoConstantSpaceBP> csbp;
  csbp = cv::cuda::createStereoConstantSpaceBP(ndisp, iters, levels, nr_plane);
  csbp->compute(d_left, d_right, d_disp);
  d_disp.download(img_disp);
}
#endif

#endif //__STEREO_COMPUTE_H__

