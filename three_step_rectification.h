#ifndef __THREE_STEP_RECTIFICATION_H__
#define __THREE_STEP_RECTIFICATION_H__

#include "mio/altro/casting.h"
#include "mio/altro/opencv.h"
#include "mio/altro/error.h"


typedef struct SThreeStepRectData{
  cv::Mat H1[2], H2[2], H3[2], F, E;
  cv::Mat R1[2], R2[2], R3[2];
  cv::Mat P[2];
  cv::Mat R[2]; //R3[i] * R2[i] * R1[i]
} threeStepRectData_t;


//return an Antisymmetric (Skew-symmetric) Matrix in 3D
void SkewSymmetric(const cv::Mat &m, cv::Mat &out){
  STD_INVALID_ARG_E( (m.rows == 1 && m.cols != 3) || (m.rows == 3 && m.cols == 1) )
  const int depth = m.depth();
  STD_INVALID_ARG_E(depth == CV_32F || depth == CV_64F)
#define SKEW_SYMMETRIC(data_type)\
    out.create(3, 3, depth);\
    data_type *md = mio::StaticCastPtr<data_type>(m.data),\
              *od = mio::StaticCastPtr<data_type>(out.data);\
    od[0] =   0.0f; od[1] = -md[2]; od[2] =  md[1];\
    od[3] =  md[2]; od[4] =   0.0f; od[5] = -md[0];\
    od[6] = -md[1]; od[7] =  md[0]; od[8] =   0.0f;
  if(depth == CV_32F){
    SKEW_SYMMETRIC(float)
  }
  else{
    SKEW_SYMMETRIC(double)
  }
#undef SKEW_SYMMETRIC
}


//both left and right epipoles have the form: epipole = [e_x, e_y, 1]^T
void GetEpipoles(const cv::Mat &F, cv::Mat &left_epipole, cv::Mat &right_epipole){
  cv::Mat A, x, b, temp;
  //F * left_epipole = 0
  A = F.colRange(0, 2);
  b = F.col(2) * -1.0;
  cv::solve(A, b, x, cv::DECOMP_SVD);
  left_epipole = cv::Mat::ones( 3, 1, x.depth() );
  temp = left_epipole.rowRange(0, 2);
  x.copyTo(temp);

  //right_epipole * F.t() = 0
  cv::Mat F_t = F.t();
  A = F_t.colRange(0, 2);
  b = F_t.col(2) * -1.0;
  cv::solve(A, b, x, cv::DECOMP_SVD);
  right_epipole = cv::Mat::ones( 3, 1, x.depth() );
  temp = right_epipole.rowRange(0, 2);
  x.copyTo(temp);
}


//compute [x1, y1, w] = H * [x0, y0, 1], then x2 = x1/w, y2 = y1/w
#define APPLY_HOMOGRAPHY(data_type, h_data, x_in, y_in, x_out, y_out)\
  {\
    const data_type w_inv = 1. / (h_data[6]*x_in + h_data[7]*y_in + h_data[8]);\
    x_out = (h_data[0]*x_in + h_data[1]*y_in + h_data[2]) * w_inv;\
    y_out = (h_data[3]*x_in + h_data[4]*y_in + h_data[5]) * w_inv;\
  }


template <typename PNT_T>
inline PNT_T ApplyHomography(const cv::Mat &H, const PNT_T &pnt_in, PNT_T &pnt_out){
  const int depth = H.depth();
  STD_INVALID_ARG_E(depth == CV_32F || depth == CV_64F)
  if(depth == CV_32F){
    const float *H_data = mio::StaticCastPtr<float>(H.data);
    APPLY_HOMOGRAPHY(float, H_data, pnt_in.x, pnt_in.y, pnt_out.x, pnt_out.y)
  }
  else{
    const double *H_data = mio::StaticCastPtr<double>(H.data);
    APPLY_HOMOGRAPHY(double, H_data, pnt_in.x, pnt_in.y, pnt_out.x, pnt_out.y)
  }
}


void ImgCentering(const cv::Size img_size, const cv::Mat &H0_src, const cv::Mat &H1_src,
                  cv::Mat &H0_dst, cv::Mat &H1_dst){
  STD_INVALID_ARG_E(H0_src.rows == 3 && H0_src.cols == 3 && H1_src.rows == 3 && H1_src.cols == 3)
  STD_INVALID_ARG_E( H0_src.type() == H1_src.type() )
  const int depth = H0_src.depth();
  STD_INVALID_ARG_E(depth == CV_32F || depth == CV_64F)
  cv::Mat extra_T0, extra_T1;
#define IMG_CENTERING(data_type, pnt_type)\
  const data_type w = static_cast<data_type>(img_size.width) * 0.5,\
                  h = static_cast<data_type>(img_size.height) * 0.5;\
  pnt_type pix[2], pix_delta[2];\
  const data_type *Hl_data = mio::StaticCastPtr<data_type>(H0_src.data),\
                  *Hr_data = mio::StaticCastPtr<data_type>(H1_src.data);\
  APPLY_HOMOGRAPHY(data_type, Hl_data, w, h, pix[0].x, pix[0].y)\
  APPLY_HOMOGRAPHY(data_type, Hr_data, w, h, pix[1].x, pix[1].y)\
  pix_delta[0].x = w - pix[0].x;\
  pix_delta[1].x = w - pix[1].x;\
  pix_delta[0].y = pix_delta[1].y = h - pix[0].y;\
  extra_T0 = (cv::Mat_<data_type>(3, 3) << 1, 0, pix_delta[0].x,\
                                           0, 1, pix_delta[0].y,\
                                           0, 0, 1);\
  extra_T1 = (cv::Mat_<data_type>(3, 3) << 1, 0, pix_delta[1].x,\
                                           0, 1, pix_delta[1].y,\
                                           0, 0, 1);
  if(depth == CV_32F){
    IMG_CENTERING(float, cv::Point2f)
  }
  else{
    IMG_CENTERING(double, cv::Point2d)
  }
#undef IMG_CENTERING
  H0_dst = extra_T0 * H0_src;
  H1_dst = extra_T1 * H1_src;
}


//applies a homography to an image
void imwarp_(const cv::Mat &src_img, cv::Mat &dst_img, const cv::Mat &H){
  STD_INVALID_ARG_E(H.rows == 3 && H.cols == 3)
  const int depth = H.depth();
  STD_INVALID_ARG_E(depth == CV_32F || depth == CV_64F)
  const cv::Mat H_inv = H.inv(); //notice we are taking the inverse here
  cv::Mat map_x(src_img.size(), CV_32F), map_y(src_img.size(), CV_32F);
  float *map_x_data = mio::StaticCastPtr<float>(map_x.data);
  float *map_y_data = mio::StaticCastPtr<float>(map_y.data);
  const size_t num_cols = src_img.cols, num_rows = src_img.rows;
  size_t idx = 0;
  //below we , where [x2, y2] (dst) is the transformation mapping coordinate for [x0, y0] (src)
#define IMWARP_MACRO(data_type)\
  const data_type *H_data = mio::StaticCastPtr<data_type>(H_inv.data);\
  for(size_t j = 0; j < num_rows; ++j){\
    const float y_ = static_cast<float>(j);\
    for(size_t i = 0; i < num_cols; ++i){\
      const float x_ = static_cast<float>(i);\
      APPLY_HOMOGRAPHY(data_type, H_data, x_, y_, map_x_data[idx], map_y_data[idx])\
      idx++;\
    }\
  }
  if(depth == CV_32F){
    IMWARP_MACRO(float)
  }
  else{
    IMWARP_MACRO(double)
  }
#undef IMWARP_MACRO

  cv::remap(src_img, dst_img, map_x, map_y, cv::INTER_LINEAR);
}

#undef APPLY_HOMOGRAPHY


void ThreeStepRectR(const cv::Mat &a, const cv::Mat &b, cv::Mat &R){
  cv::Mat t, t_ss;
  double cos_theta = a.dot(b) / ( cv::norm(a) * cv::norm(b) );
  cv::normalize(a.cross(b), t);
  if(cos_theta < 0){
    cos_theta *= -1.0;
    t *= -1.0;
  }
  SkewSymmetric(t, t_ss);
  const cv::Mat I = cv::Mat::eye( 3, 3, a.depth() );
  const double sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
  R = I + sin_theta * t_ss + (1.0 - cos_theta) * t_ss * t_ss;
  if(false){
    printf("cos_theta: %f\n", cos_theta);
    mio::Print(t, "t");
    mio::Print(t_ss, "t_ss");
    mio::Print(R, "R");
    std::cout << "\n\n";
  }
}


void GetH1H2(const cv::Mat &K, const cv::Mat &epipole, cv::Mat &H1, cv::Mat &H2, cv::Mat &R1, cv::Mat &R2){
  const int depth = epipole.depth();
  STD_INVALID_ARG_E(depth == CV_32F || depth == CV_64F)
  const cv::Mat K_inv = K.inv();
  cv::Mat a, b, t, t_ss;

  //step 1 - the first step sends both epipoles to infinity [e_x; e_y; 0]
  cv::Mat epipole_inf = cv::Mat::zeros(3, 1, depth);
  cv::Mat epipole_temp = epipole.rowRange(0, 2), epipole_inf_temp = epipole_inf.rowRange(0, 2);
  epipole_temp.copyTo(epipole_inf_temp);
  a = K_inv * epipole;
  b = K_inv * epipole_inf;
  ThreeStepRectR(a, b, R1);
  H1 = K * R1 * K_inv;
  if(false){
    mio::Print(epipole, "epipole");
    mio::Print(epipole_inf, "epipole_inf");
    mio::Print(H1, "H1");
  }

  //step 2 - then the epipoles are sent to infinity in the x-direction
  cv::Mat epipole_i;
  if(depth == CV_32F)
    epipole_i = (cv::Mat_<float>(3, 1) << 1, 0, 0);
  else
    epipole_i = (cv::Mat_<double>(3, 1) << 1, 0, 0);
  a = K_inv * epipole_inf;
  b = K_inv * epipole_i;
  ThreeStepRectR(a, b, R2);
  H2 = K * R2 * R1 * K_inv;
  if(false){
    mio::Print(epipole_i, "epipole_i");
    mio::Print(H2, "H2");
  }
}


//ThreeStepRectification implementation
void ResidualRotation(const cv::Mat &E, cv::Mat &R){
  STD_INVALID_ARG_E(E.depth() == CV_64F)

  cv::Mat W, U, Vt;
  cv::SVD::compute(E, W, U, Vt);

  const cv::Mat W_ = (cv::Mat_<double>(3, 3) << 1, 0, 0,
                                                0, 1, 0,
                                                0, 0, 0);
  const cv::Mat E_ = U * W_ * Vt;
  //const cv::Mat F2_ = K_r_inv.t() * E_ * K_l_inv; //not used

  STD_INVALID_ARG_E(E_.depth() == CV_64F)
  const cv::Mat det_mat = ( cv::Mat_<double>(2, 2) << E_.at<double>(2,1),  E_.at<double>(2,2),
                                                     -E_.at<double>(1,1), -E_.at<double>(1,2) );
  const double sub_det = cv::determinant(det_mat);

  R = cv::Mat::zeros( 3, 3, E.type() );
  cv::Mat R_sub( R, cv::Rect(1, 1, 2, 2) );
  det_mat.copyTo(R_sub);
  R.at<double>(0,0) = (std::acos( E_.at<double>(2,1) ) > M_PI*0.5f) ? -std::sqrt(sub_det) : std::sqrt(sub_det);
  cv::Mat t = U.col(2) * 1.0;
}


void DecomposeEssentialMat(const cv::Mat &E, cv::Mat &R1, cv::Mat &R2, cv::Mat &t){
  STD_INVALID_ARG_E(E.depth() == CV_64F)
  STD_INVALID_ARG_E(E.cols == 3 && E.rows == 3)

  cv::Mat D, U, Vt;
  cv::SVD::compute(E, D, U, Vt);
  cv::Mat W = (cv::Mat_<double>(3, 3) << 0, -1, 0,
                                         1,  0, 0,
                                         0,  0, 1);
  W.convertTo( W, E.type() );
  //two possible rotations: W or W.t()
  R1 = U * W * Vt;
  R2 = U * W.t() * Vt;
  t = U.col(2) * 1.0;
}


void ThreeStepRect(const stereoCalData_t &cal, threeStepRectData_t &rect,
                   const std::string rotation_type = "Standard"){
  cv::Mat epipole_l, epipole_r;
  GetEpipoles(cal.F, epipole_l, epipole_r);

  cv::Mat H1_[2], H2_[2], H3_[2];
  GetH1H2(cal.K[0], epipole_l, H1_[0], H2_[0], rect.R1[0], rect.R2[0]);
  GetH1H2(cal.K[1], epipole_r, H1_[1], H2_[1], rect.R1[1], rect.R2[1]);

  //step 3 - compute H3 - the residual rotation between both cameras around their baseline
  //is compensated to achieve the rectification
  //1. compute new F using R1 and R2
  //2. compute new E from F
  //3. extract residual rotation (R3) from new E
  cv::Mat K_inv[2];
  K_inv[0] = cal.K[0].inv();
  K_inv[1] = cal.K[1].inv();
  rect.F = K_inv[1].t() * rect.R2[1] * rect.R1[1] * cal.K[1].t() * cal.F * cal.K[0] * rect.R1[0].inv() * rect.R2[0].inv() * K_inv[0];
  rect.E = cal.K[1].t() * rect.F * cal.K[0];

  cv::Mat R1, R2, t;
  DecomposeEssentialMat(rect.E, R1, R2, t);
  if(rotation_type == "Standard")
    ResidualRotation(rect.E, rect.R3[0]);
  else if(rotation_type == "OpenCV R1")
    rect.R3[0] = R1;
  else if(rotation_type == "OpenCV R2")
    rect.R3[0] = R2;
  else
    STD_INVALID_ARG_M("Invalid rotation_type string")

  rect.R3[1] = cv::Mat::eye( 3, 3, rect.R3[0].depth() );
  H3_[0] = cal.K[0] * rect.R3[0] * rect.R2[0] * rect.R1[0] * K_inv[0];
  H3_[1] = cal.K[1] * rect.R3[1] * rect.R2[1] * rect.R1[1] * K_inv[1];
  rect.R[0] = rect.R3[0] * rect.R2[0] * rect.R1[0];
  rect.R[1] = rect.R3[1] * rect.R2[1] * rect.R1[1];

  ImgCentering(cal.img_size, H1_[0], H1_[1], rect.H1[0], rect.H1[1]);
  ImgCentering(cal.img_size, H2_[0], H2_[1], rect.H2[0], rect.H2[1]);
  ImgCentering(cal.img_size, H3_[0], H3_[1], rect.H3[0], rect.H3[1]);
}

#endif //__THREE_STEP_RECTIFICATION_H__

