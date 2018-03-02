#ifndef __STEREO_TYPES__
#define __STEREO_TYPES__

#include <array>

#if ICV_OPENCV_VERSION_MAJOR < 3
#include "opencv2/core/core.hpp"
#else
#include "opencv2/core/core.hpp"
#endif

typedef struct SStereoCalData{
  cv::Mat R, T; //extrinsic
  cv::Mat K[2], D[2]; //intrinsic
  cv::Mat E, F; //essential, fundamental
  cv::Mat H1, H2,
          R1, R2, //3x3 rectification transform (rotation matrix) //TODO - change to R1_rect or rectif_transform[]
          P1, P2, //3x4 projection matrix in the new (rectified) coordinate systems //TODO - change to P1_rect or ...
          Q; //4x4 disparity-to-depth mapping matrix ( see reprojectImageTo3D() )
  cv::Rect  valid_roi[2]; //valid image region with rectified images
  cv::Size img_size; //size of original raw image from single camera
  std::array<std::vector<std::string>, 2> good_img_file_names;
  //img_cal_pnts[first/second camera][camera image idx][pnt idx]
  std::array<std::vector< std::vector<cv::Point2f> >, 2> img_cal_pnts;
  int num_img_per_cam;

  void SetNumValidImgPair(const int num){
    num_img_per_cam = num;
  }

  int GetNumImgPerCam() const{
    return num_img_per_cam;
  }

  void ClearRectData(){
    H1.release();
    H2.release();
    R1.release();
    R2.release();
    P1.release();
    P2.release();
    Q.release();
    valid_roi[0] = cv::Rect();
    valid_roi[1] = cv::Rect();
  }

  void Print(const bool single_cam = false){
    std::cout << "\nIntrinsic Parameters\n\n";
    if(single_cam){
      std::cout << "M1:\n" << K[0] << "\n\n"
                   "D1:\n" << D[0] << "\n";
    }
    else{
      std::cout << "M1:\n" << K[0] << "\n\n"
                   "M2:\n" << K[1] << "\n\n"
                   "D1:\n" << D[0] << "\n\n"
                   "D2:\n" << D[1] << "\n";

      std::cout << "\nExtrinsic Parameters\n\n" << 
                   "R:\n" << R << "\n\n"
                   "T:\n" << T << "\n";

      std::cout << "\nFundamental/Essential Matrices\n\n" << 
                   "F:\n" << F << "\n\n"
                   "E:\n" << E << "\n";

      std::cout << "\nRectification Matrices:\n\n" <<
                   "H1:\n" << H1 << "\n\n" <<
                   "H2:\n" << H2 << "\n\n" <<
                   "R1:\n" << R1 << "\n\n" <<
                   "R2:\n" << R2 << "\n\n" <<
                   "P1:\n" << P1 << "\n\n" <<
                   "P2:\n" << P2 << "\n\n" <<
                   "Q:\n"  << Q  << "\n\n" <<
                   "valid_roi[0]:\n" << valid_roi[0] << "\n\n" <<
                   "valid_roi[1]:\n" << valid_roi[1] << "\n";
    }

    std::cout << "\nnum_img_per_cam: " << num_img_per_cam << "\n";

    std::cout << "\ncamera_0 file names:\n";
    for(auto &str : good_img_file_names[0])
      std::cout << str << " ";
    if(!single_cam){
      std::cout << "\n\ncamera_1 file names:\n";
      for(auto &str : good_img_file_names[1])
        std::cout << str << " ";
   }
    std::cout << "\n\n";

    std::cout << "img_cal_pnts[0].size(): " << img_cal_pnts[0].size() << " - ";
    for(auto &vec : img_cal_pnts[0])
      std::cout << vec.size() << " ";
    if(!single_cam){
      std::cout << "\n\nimg_cal_pnts[1].size(): " << img_cal_pnts[1].size() << " - ";
      for(auto &vec : img_cal_pnts[1])
        std::cout << vec.size() << " ";
    }
    std::cout << "\n\n";
  }

} stereoCalData_t;


typedef struct SCameraCalTarget{
  std::string type_str;
  cv::Size size;
  float spacing;
  std::vector<cv::Point3f> point_vec;

  SCameraCalTarget(){}

  SCameraCalTarget(std::string type_str_, cv::Size size_, float spacing_) : 
                   type_str(type_str_), size(size_), spacing(spacing_) {
    std::cout << "Target type: " << type_str << std::endl;
    printf("target dimensions (height, width): (%d, %d)\n", size.height, size.width);
    printf("targetSpacing: %f\n", spacing);

    //create the virtual calibration target
    point_vec.resize(size.height*size.width);
    if(type_str == "a-circle"){
      const double SQRT3 = 1.73205080757;
      cv::Point2f asymm_target_spacing(spacing/2.0*SQRT3, spacing/2.0);
      for(size_t i = 0, r = 0; r < size.height; r++)
        for(size_t c = 0; c < size.width; c++){
          //equally spaced point asymmetric pattern
          //point_vec[i++] = cv::Point3f( c*asymm_target_spacing.x, (2*r + c%2)*asymm_target_spacing.x, 0 );

          //taken from modules/calib3d/src/circlesgrid.cpp
          //point_vec[i++] = cv::Point3f( (2*c + r%2)*target.spacing, r*target.spacing, 0 );

          //taken from modules/calib3d/src/circlesgrid.cpp, equally spaced
          point_vec[i++] = cv::Point3f( (2*c + r%2)*asymm_target_spacing.x, r*asymm_target_spacing.y, 0 );
        }
    }
    else
      for(size_t i = 0, r = 0; r < size.height; r++)
        for(size_t c = 0; c < size.width; c++)
          point_vec[i++] = cv::Point3f(c*spacing, r*spacing, 0);
  }
} camCalTarget_t;

#endif //__STEREO_TYPES__

