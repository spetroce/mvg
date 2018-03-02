#ifndef __STEREO_IO__
#define __STEREO_IO__

#include <QFile>
#include <QXmlStreamWriter>
#include <QMessageBox>
#include <iostream>
#include <unordered_map>

#include "opencv2/core.hpp"

#include "mio/qt/qt_xml.h"
#include "stereo_types.h"
#include "mio/altro/opencv.h"


inline void SaveCalibrationPoints(const std::string file_path, const stereoCalData_t &cal_data){
  QFile file(QString::fromStdString(file_path) + "/img_cal_pnts.xml");
  if( !file.open(QIODevice::WriteOnly) ) //open a file
    QMessageBox::warning(0, "Read only", "The file is in read only mode");
  else{
    QXmlStreamWriter xml;
    xml.setDevice(&file); //set device (here file)to streamwriter
    xml.writeStartDocument(); //write XML version number

    xml.writeStartElement("img_cal_pnts");
    xml.writeTextElement( "num_valid_img_pair", INT_TO_QSTR( cal_data.GetNumImgPerCam() ) );

    #define WRITE_IMG_CAL_PNTS_XML(camera_name_str, camera_index)\
    xml.writeStartElement(camera_name_str);\
    for(size_t i = 0; i < cal_data.GetNumImgPerCam(); ++i){\
      xml.writeStartElement("cal_image");\
      const QString file_name = QString::fromStdString(cal_data.good_img_file_names[camera_index][i]);\
      xml.writeAttribute( "file_name",  file_name);\
      for(auto &cv_point2f : cal_data.img_cal_pnts[camera_index][i]){\
        xml.writeStartElement("point2f");\
        xml.writeAttribute( "x", QString().sprintf("%.6f", cv_point2f.x) );\
        xml.writeAttribute( "y", QString().sprintf("%.6f", cv_point2f.y) );\
        xml.writeEndElement(); /*"point2f"*/\
      }\
      xml.writeEndElement(); /*"calibration_image"*/\
    }\
    xml.writeEndElement(); /*camera_name_str*/

    WRITE_IMG_CAL_PNTS_XML("camera_0", 0)
    WRITE_IMG_CAL_PNTS_XML("camera_1", 1)
    #undef WRITE_IMG_CAL_PNTS_XML

    xml.writeEndElement(); //"img_cal_pnts"
  }
  file.close();
}


inline void SaveCameraCalData(const std::string file_path, const std::string intrinsic_file_name,
                              const stereoCalData_t &cal_data){
  try{
    cv::FileStorage fs;
    //save intrinsic parameters
    fs.open(file_path + "/" + intrinsic_file_name + ".yml", cv::FileStorage::WRITE);
    if( fs.isOpened() ){
      fs << "M" << cal_data.K[0] << "D" << cal_data.D[0] << "ImageSize" << cal_data.img_size;
      fs.release();
    }
    else
      std::cout << "Error: can not save the intrinsic parameters\n";

    SaveCalibrationPoints(file_path, cal_data);
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }
}


inline void SaveStereoCalData(const std::string file_path, const std::string intrinsic_file_name,
                              const std::string extrinsic_file_name, const stereoCalData_t &cal_data,
                              const bool only_rectification_params = false){
  cv::FileStorage fs;

  if(!only_rectification_params){
    try{
      //save intrinsic parameters
      fs.open(file_path + "/" + intrinsic_file_name + ".yml", cv::FileStorage::WRITE);
      if( fs.isOpened() ){
        fs << "M1" << cal_data.K[0] << "D1" << cal_data.D[0] << 
              "M2" << cal_data.K[1] << "D2" << cal_data.D[1] << "ImageSize" << cal_data.img_size;
        fs.release();
      }
      else
        std::cout << "Error: can not save the intrinsic parameters\n";
    }
    catch(cv::Exception &e){
      FRIENDLY_RETHROW(e)
    }

    try{
      //save all extrinsic parameters
      fs.open(file_path + "/" + extrinsic_file_name + ".yml", cv::FileStorage::WRITE);
      if( fs.isOpened() ){
        fs << "R" << cal_data.R << "T" << cal_data.T;
        fs.release();
      }
      else
        std::cout << "Error: can not save the extrinsic parameters\n";
    }
    catch(cv::Exception &e){
      FRIENDLY_RETHROW(e)
    }

    try{
      //save fundamental/essential matrices
      fs.open(file_path + "/" + std::string("fundamentalEssential.yml"), cv::FileStorage::WRITE);
      if( fs.isOpened() ){
        fs << "F" << cal_data.F << "E" << cal_data.E;
        fs.release();
      }
      else
        std::cout << "Error: can not save the fundamental/essential matrices\n";
    }
    catch(cv::Exception &e){
      FRIENDLY_RETHROW(e)
    }

    SaveCalibrationPoints(file_path, cal_data);
  }

  try{
    //save rectification parameters
    fs.open(file_path + "/" + extrinsic_file_name + "_rect.yml", cv::FileStorage::WRITE);
    if( fs.isOpened() ){
      fs << "R1" << cal_data.R1 << "R2" << cal_data.R2 <<
            "P1" << cal_data.P1 << "P2" << cal_data.P2 << "Q" << cal_data.Q <<
            "validPixROI1" << cal_data.valid_roi[0] << "validPixROI2" << cal_data.valid_roi[1];
      fs.release();
    }
    else
      std::cout << "Error: can not save the rectified extrinsic parameters\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }
}


inline cv::Point2f ParseXMLPoint2f(const QXmlStreamReader &xml){
  QXmlStreamAttributes pnt_attr = xml.attributes();
  float x, y;
  if(pnt_attr.hasAttribute("x") )
    x = pnt_attr.value("x").toString().toFloat();
  if(pnt_attr.hasAttribute("y") )
    y = pnt_attr.value("y").toString().toFloat();
  return cv::Point2f(x, y);
}


inline void ParseXMLCamera(QXmlStreamReader &xml, const QString &camera_name,
                           const size_t cam_idx, stereoCalData_t &cal_data){
  if(xml.name() == camera_name){
    while( !xml.atEnd() && !xml.hasError() ){
      xml.readNext();
      if(xml.isEndElement() && xml.name() == camera_name)
        break;
      if(xml.isStartElement() && xml.name() == "cal_image"){
        QXmlStreamAttributes attr = xml.attributes();
        if(attr.hasAttribute("file_name") )
          cal_data.good_img_file_names[cam_idx].push_back( attr.value("file_name").toString().toStdString() );
        std::vector<cv::Point2f> point_vec;
        while( !xml.atEnd() && !xml.hasError() ){
          xml.readNext();
          if(xml.isEndElement() && xml.name() == "cal_image")
            break;
          if(xml.isStartElement() && xml.name() == "point2f")
            point_vec.push_back( ParseXMLPoint2f(xml) );
        }
        cal_data.img_cal_pnts[cam_idx].push_back(point_vec);
      }
    }
  }
}


inline void LoadCameraCalData(const std::string file_path, const std::string intrinsic_file_name_1,
                              const std::string intrinsic_file_name_2, stereoCalData_t &cal_data){
  cv::FileStorage fs;
  std::string file_full;
  try{
    //load intrinsic parameters for first camera
    file_full = file_path + "/" + intrinsic_file_name_1 + ".yml";
    fs.open(file_full, cv::FileStorage::READ);
    if( fs.isOpened() ){
      fs["M"] >> cal_data.K[0];
      fs["D"] >> cal_data.D[0];
      fs["ImageSize"] >> cal_data.img_size;
      fs.release();
    }
    else
      std::cout << "Error: can not load the intrinsic parameters from file " << file_full << "\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }

  try{
    //load intrinsic parameters for second camera
    file_full = file_path + "/" + intrinsic_file_name_2 + ".yml";
    fs.open(file_full, cv::FileStorage::READ);
    if( fs.isOpened() ){
      fs["M"] >> cal_data.K[1];
      fs["D"] >> cal_data.D[1];
      fs["ImageSize"] >> cal_data.img_size;
      fs.release();
    }
    else
      std::cout << "Error: can not load the intrinsic parameters from file " << file_full << "\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }
}


inline void LoadStereoCalData(const std::string file_path, const std::string intrinsic_file_name,
                              const std::string extrinsic_file_name, stereoCalData_t &cal_data/*,
                              std::vector<cv::Mat> &cal_img_vec*/){
  cv::FileStorage fs;
  try{
    //save intrinsic parameters
    fs.open(file_path + "/" + intrinsic_file_name + ".yml", cv::FileStorage::READ);
    if( fs.isOpened() ){
      fs["M1"] >> cal_data.K[0];
      fs["D1"] >> cal_data.D[0];
      fs["M2"] >> cal_data.K[1];
      fs["D2"] >> cal_data.D[1];
      fs["ImageSize"] >> cal_data.img_size;
      fs.release();
    }
    else
      std::cout << "Error: can not load the intrinsic parameters\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }

  try{
    //save all extrinsic parameters
    fs.open(file_path + "/" + extrinsic_file_name + ".yml", cv::FileStorage::READ);
    if( fs.isOpened() ){
      fs["R"] >> cal_data.R;
      fs["T"] >> cal_data.T;
      fs.release();
    }
    else
      std::cout << "Error: can not load the extrinsic parameters\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }

  try{
    //save fundamental/essential matrices
    fs.open(file_path + "/" + std::string("fundamentalEssential.yml"), cv::FileStorage::READ);
    if( fs.isOpened() ){
      fs["F"] >> cal_data.F;
      fs["E"] >> cal_data.E;
      fs.release();
    }
    else
      std::cout << "Error: can not load the fundamental matrix\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }

  try{
    //save rectification parameters
    fs.open(file_path + "/" + extrinsic_file_name + "_rect.yml", cv::FileStorage::READ);
    if( fs.isOpened() ){
      fs["R1"] >> cal_data.R1;
      fs["R2"] >> cal_data.R2;
      fs["P1"] >> cal_data.P1;
      fs["P2"] >> cal_data.P2;
      fs["Q"] >> cal_data.Q;
      fs["validPixROI1"] >> cal_data.valid_roi[0];
      fs["validPixROI2"] >> cal_data.valid_roi[1];
      fs.release();
    }
    else
      std::cout << "Error: can not load the rectified extrinsic parameters\n";
  }
  catch(cv::Exception &e){
    FRIENDLY_RETHROW(e)
  }

  cal_data.good_img_file_names[0].clear();
  cal_data.good_img_file_names[1].clear();
  QFile file(QString::fromStdString(file_path) + "/img_cal_pnts.xml");
  if( !file.open(QIODevice::ReadOnly | QIODevice::Text) )
    QMessageBox::warning(0, "LoadStereoCalData() Error", "Could not open file img_cal_pnts.xml");
  else{
    QXmlStreamReader xml(&file);
    while( !xml.atEnd() && !xml.hasError() ){
      QXmlStreamReader::TokenType token = xml.readNext();
      if( xml.isStartDocument() )
        continue;
      if( xml.isStartElement() ){
        if(xml.name() == "img_cal_pnts")
          continue;
        if(xml.name() == "num_valid_img_pair"){
          cal_data.SetNumValidImgPair( xml.readElementText().toInt() );
          cal_data.img_cal_pnts[0].reserve( cal_data.GetNumImgPerCam() );
          cal_data.img_cal_pnts[1].reserve( cal_data.GetNumImgPerCam() );       
        }
        ParseXMLCamera(xml, "camera_0", 0, cal_data);
        ParseXMLCamera(xml, "camera_1", 1, cal_data);
      }
    }
  }
}


inline cv::Mat LoadCalImgDefaultProc(const std::string file_full, const size_t bit_depth = 8){
  STD_INVALID_ARG_E(bit_depth == 8 || bit_depth == 10 || bit_depth == 12 || bit_depth == 14 || bit_depth == 16)
  cv::Mat src = mio::ReadImage(file_full, cv::IMREAD_GRAYSCALE), dst;
  if(bit_depth == 8)
    return src;
  else{
    std::unordered_map<size_t, float> scale_down_map { {10, 0.249f},        // (2^8 - 1) ÷ (2^10 - 1)
                                                       {12, 0.0622f},       // (2^8 - 1) ÷ (2^12 - 1)
                                                       {14, 0.01556491f},   // (2^8 - 1) ÷ (2^14 - 1)
                                                       {16, 0.00389105f} }; // (2^8 - 1) ÷ (2^16 - 1)
    src.convertTo(dst, CV_8U, scale_down_map[bit_depth]);
    return dst;
  }
  return dst;
}


inline void LoadCalImages(const std::string &file_path, const camCalTarget_t &cal_target,
                          const stereoCalData_t &cal_data, std::vector<cv::Mat> &cal_img_vec, const size_t num_camera){
  STD_INVALID_ARG_E(num_camera == 1 || num_camera == 2)
  cal_img_vec.resize( cal_data.GetNumImgPerCam() );
  cv::Mat color_cal_img[2];

  for(size_t i = 0; i < cal_data.GetNumImgPerCam(); ++i){
    for(size_t k = 0; k < num_camera; ++k){
      //draw detected calibration target features on cal image
      cv::Mat img = LoadCalImgDefaultProc(file_path + "/" + cal_data.good_img_file_names[k][i]), color_img;
      cv::cvtColor(img, color_img, cv::COLOR_GRAY2BGR);
      cv::drawChessboardCorners(color_img, cal_target.size, cal_data.img_cal_pnts[k][i], true);
      double sf = 480./MAX(img.rows, img.cols);
      cv::resize(color_img, color_cal_img[k], cv::Size(), sf, sf);
    }

    //save a combined left/right cal image for viewing
    if(num_camera == 2){
      const size_t rows = color_cal_img[0].rows, cols = color_cal_img[0].cols;
      cv::Mat dst( rows, cols*2, color_cal_img[0].type() );
      cv::Mat left = dst( cv::Rect(0, 0, cols, rows) ),
              right = dst( cv::Rect(cols, 0, cols, rows) );
      color_cal_img[0].copyTo(left);
      color_cal_img[1].copyTo(right);
      cal_img_vec[i] = dst.clone();
    }
    else
      cal_img_vec[i] = color_cal_img[0].clone();
  }
}

#endif //__STEREO_IO__
