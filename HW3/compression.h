/*
 *	2015/11/05
 *	Tony Guo
 *	R04458006
 *
 *	HW3 - DCT IDCT Quantization InverseQuantization RLD RLC HuffmandVLC HuffmanVLD ArithmeticVLC ArithmeticVLD
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

#include "ac.h"

struct HuffmanNode
{
   HuffmanNode(int prob,char data,std::string code,HuffmanNode *L,HuffmanNode *R):
      prob(prob),data(data),code(code),L(L),R(R){}
   int prob;
   char data;
   std::string code;
   HuffmanNode *L;
   HuffmanNode *R;
   bool operator==(const char data)const
   {
      return (data == this->data);
   }
};

inline bool HuffmanArraySort(HuffmanNode* &l,HuffmanNode* &r){return l->prob < r->prob;}

class Compression
{
   public:
      Compression(std::string path,unsigned int raw_width,unsigned int raw_height,unsigned int block_size = 8):
	 path(path),raw_width(raw_width),raw_height(raw_height),block_size(block_size){LoadRawImage();}
      void LoadRawImage();
      void OneD_Block_DCT();
      void OneD_Block_IDCT();
      void Quantization();
      void IQuantization();
      void RLC();
      void RLD();
      void HuffmanVLC();
      void HuffmanVLD();
      void AriVLC();
      void AriVLD();
      double GetHuffmanByte(){return huffman_total_byte;}
      double GetAriByte(){return ari_total_byte;}

      double PSNRComputing()
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

      double C(double value)
      {
	 return (value == 0)?1.0f/std::sqrt(2.0f):1.0f;
      }

      double limit(double value)
      {
	 return (value < 0)?0:(value > 255)?255:value;
      }
   private:
      void CreateHuffmanTree(std::map<char,int>&,std::vector<HuffmanNode*>&);
      void CreateHuffmanCode(std::vector<HuffmanNode>&,HuffmanNode*,std::string);

      std::string path;
      unsigned int raw_width,raw_height;
      unsigned int block_size;
      std::vector<std::vector<unsigned char> > origin_image;
      std::vector<std::vector<double> > dct;
      std::vector<std::vector<unsigned char> > idct;


      double huffman_total_byte;
      double ari_total_byte;
};

#endif
