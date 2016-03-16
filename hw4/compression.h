/*
 *	2015/11/22
 *	Tony Guo
 *	R04458006
 *
 *	HW4
 *
 */

#ifndef _COMPRESSION_
#define _COMPRESSION_

#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <string>
#include <vector>
#include <bitset>
#include <map>

#include <cmath>

#include <memory>


enum{
   FULL = 0,
   FAST
};

enum{
   INTRA = 0,
   INTER
};

class Compression
{
   public:
      Compression(std::string path,unsigned int raw_width,unsigned int raw_height,unsigned int block_size = 8):
	 path(path),raw_width(raw_width),raw_height(raw_height),block_size(block_size){}
      void LoadRawImage(std::vector<std::vector<unsigned char> >&);
      void LoadRawVideo(std::vector<std::vector<std::vector<unsigned char> > >&);
      void OneD_Block_DCT(std::vector<std::vector<unsigned char> >&,std::vector<std::vector<double> >&);
      void OneD_Block_IDCT(std::vector<std::vector<double> >&,std::vector<std::vector<unsigned char> >&);
      void Quantization(std::vector<std::vector<double> >&,const int);
      void IQuantization(std::vector<std::vector<double> >&,const int);
      void NextFrame(std::vector<std::vector<std::vector<unsigned char> > >&,std::vector<std::vector<unsigned char> >&,unsigned int);
      void ME(std::vector<std::vector<unsigned char> >&,std::vector<std::vector<unsigned char> >&,int,unsigned int);
      void MC(std::vector<std::vector<unsigned char> >&,std::vector<std::vector<unsigned char> >&);

      void SubImage(std::vector<std::vector<unsigned char> > &src1,std::vector<std::vector<unsigned char> > &src2,std::vector<std::vector<unsigned char> > &dst)
      {
	 for(size_t i = 0;i < src1.size();++i){
	    for(size_t j = 0;j < src1[i].size();++j){
	       dst[i][j] = limit(src1[i][j] - src2[i][j]);
	    }
	 }
      }
      void AddImage(std::vector<std::vector<unsigned char> > &src1,std::vector<std::vector<unsigned char> > src2,std::vector<std::vector<unsigned char> > dst)
      {
	 for(size_t i = 0;i < src1.size();++i){
	    for(size_t j = 0;j < src1[i].size();++j){
	       dst[i][j] = limit(src1[i][j] + src2[i][j]);
	    }
	 }
      }

      inline unsigned int TotalFrame(){return total_frame;}

      inline double PSNRComputing(std::vector<std::vector<unsigned char> > &origin_image,std::vector<std::vector<unsigned char> > &idct)
      {
	 double mse = 0,psnr = 0;

	 for(int i = 0;i < raw_width;i++){
	    for(int j = 0;j < raw_height;j++){
	       mse += (origin_image[i][j] - idct[i][j]) * (origin_image[i][j] - idct[i][j]);
	    }
	 }

	 if(mse != 0){
	    mse /= raw_width * raw_height;
	    psnr = 10 * std::log10((255 * 255) / mse);
	 }
	 else{
	    psnr = 999999;
	 }

	 //std::cout << "PSNR : " << psnr << " dB" << std::endl;
	 return psnr;
      }

      inline double C(double value)
      {
	 return (value == 0)?1.0f/std::sqrt(2.0f):1.0f;
      }

      inline double limit(double value)
      {
	 return (value < 0)?0:(value > 255)?255:value;
      }
   private:
      std::string path;
      unsigned int raw_width,raw_height;
      unsigned int block_size;
      unsigned int total_frame;

      std::vector< std::pair<int,int> > MV;
      /*
      std::vector<std::vector<unsigned char> > origin_image;
      std::vector<std::vector<double> > dct;
      std::vector<std::vector<unsigned char> > idct;

      std::vector<std::vector<std::vector<unsigned char> > > video;
      */
};

#endif
