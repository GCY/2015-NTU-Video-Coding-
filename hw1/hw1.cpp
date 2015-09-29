/*
 *	2015/09/29
 *	Tony Guo
 *	R04458006
 *
 *
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>

const int RAW_WIDTH = 512;
const int RAW_HEIGHT = 512;
const int RAW_CHANNEL = 3;

inline double normalize(double value,double in_min,double in_max,double out_min,double out_max)
{
  return (value - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

inline double limit(double value)
{
   if(value > 255){return 255;}
   else if(value < 0){return 0;}

   return value;
}

void RGBChannelSplit(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *r,char *g,char *b)
{
   std::fstream r_channel(r,std::ios::binary | std::ios::out);
   std::fstream g_channel(g,std::ios::binary | std::ios::out);
   std::fstream b_channel(b,std::ios::binary | std::ios::out);
   
   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 r_channel << image[i][j][0];
	 g_channel << image[i][j][1];
	 b_channel << image[i][j][2];
      }
   }

   r_channel.close();
   g_channel.close();
   b_channel.close();

}

void Y(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *y)
{
   std::fstream y_channel(y,std::ios::binary | std::ios::out);

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 y_channel << (unsigned char)((image[i][j][0] * 0.299f) + (image[i][j][1] * 0.587f) + (image[i][j][2] * 0.114f));
      }
   }

   y_channel.close();
}

void I(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *i)
{
   std::fstream i_channel(i,std::ios::binary | std::ios::out);

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 i_channel << (unsigned char)(normalize(((image[i][j][0] * 0.596f) + (image[i][j][1] * -0.275f) + (image[i][j][2] * -0.321f)),-127,127,0,255) );
      }
   }

   i_channel.close();
}

void Q(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *q)
{
   std::fstream q_channel(q,std::ios::binary | std::ios::out);

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 q_channel << (unsigned char)(normalize(((image[i][j][0] * 0.212f) + (image[i][j][1] * -0.523f) + (image[i][j][2] * 0.311f)),-127,127,0,255) );
      }
   }

   q_channel.close();
}

void U(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *u)
{
   std::fstream u_channel(u,std::ios::binary | std::ios::out);

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 u_channel << (unsigned char)(normalize(((image[i][j][0] * -0.147f) + (image[i][j][1] * -0.289f) + (image[i][j][2] * 0.436f)),-111,111,0,255) );
      }
   }

   u_channel.close();
}

void V(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *v)
{
   std::fstream v_channel(v,std::ios::binary | std::ios::out);

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 v_channel << (unsigned char)(normalize(((image[i][j][0] * 0.615f) + (image[i][j][1] * -0.515f) + (image[i][j][2] * -0.100f)),-157,157,0,255) );
      }
   }

   v_channel.close();
}

void YUV420ToRGBandPSNR(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL],char *y,char *u,char *v)
{
   std::fstream y_channel(y,std::ios::binary | std::ios::in);
   std::fstream u_channel(u,std::ios::binary | std::ios::in);
   std::fstream v_channel(v,std::ios::binary | std::ios::in);
   std::fstream yuv420_channel("yub420.raw",std::ios::binary | std::ios::out);
   std::fstream rgb420_channel("rgb420.raw",std::ios::binary | std::ios::out);

   unsigned char y_image[RAW_WIDTH][RAW_HEIGHT];
   unsigned char u_image[RAW_WIDTH][RAW_HEIGHT];
   unsigned char v_image[RAW_WIDTH][RAW_HEIGHT];
   unsigned char u420_image[RAW_WIDTH/2][RAW_HEIGHT/2];
   unsigned char v420_image[RAW_WIDTH/2][RAW_HEIGHT/2];

   unsigned char rgb420_image[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL];


   unsigned char data;
   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 y_channel.read((char*) &data, sizeof(char));
	 y_image[i][j] = data;
      }
   }

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 u_channel.read((char*) &data, sizeof(char));
	 u_image[i][j] = data;
      }
   }

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 v_channel.read((char*) &data, sizeof(char));
	 v_image[i][j] = data;
      }
   }

   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 yuv420_channel << y_image[i][j];
      }
   }

   for(int i = 0,w = 0;i < RAW_WIDTH;i += 2,++w){
      for(int j = 0,h = 0;j < RAW_HEIGHT;j += 2,++h){
	 yuv420_channel << u_image[i][j];
	 u420_image[w][h] = u_image[i][j];
      }
   }

   for(int i = 0,w = 0;i < RAW_WIDTH;i += 2,++w){
      for(int j = 0,h = 0;j < RAW_HEIGHT;j += 2,++h){
	 yuv420_channel << v_image[i][j];
	 v420_image[w][h] = v_image[i][j];
      }
   }
   
   for(int i = 0;i < RAW_WIDTH;++i){
      for(int j = 0;j < RAW_HEIGHT;++j){
	 double normalize_u = normalize(u420_image[i/2][j/2],0,255,-111,111);
	 double normalize_v = normalize(v420_image[i/2][j/2],0,255,-157,157);
	 rgb420_channel << (unsigned char)limit(y_image[i][j] + normalize_v * 1.140);
	 rgb420_channel << (unsigned char)limit(y_image[i][j] + normalize_u * -0.395 + normalize_v * -0.581); 
	 rgb420_channel << (unsigned char)limit(y_image[i][j] + normalize_u * 2.032);

	 rgb420_image[i][j][0] = (unsigned char)limit(y_image[i][j] + normalize_v * 1.140);
	 rgb420_image[i][j][1] = (unsigned char)limit(y_image[i][j] + normalize_u * -0.395 + normalize_v * -0.581); 
	 rgb420_image[i][j][2] = (unsigned char)limit(y_image[i][j] + normalize_u * 2.032);

      }
   }

   double mse = 0,psnr = 0;

   for(int i = 0;i < RAW_WIDTH;i++){
      for(int j = 0;j < RAW_HEIGHT;j++){
	 for(int k = 0;k < RAW_CHANNEL;k++){
	    mse += (image[i][j][k] - rgb420_image[i][j][k]) * (image[i][j][k] - rgb420_image[i][j][k]);
	 } 
      }
   }

   mse /= RAW_WIDTH * RAW_HEIGHT * RAW_CHANNEL;
   psnr = 10 * std::log10((255 * 255) / mse);

   std::cout << "PSNR : " << psnr << " dB" << std::endl;

   y_channel.close();
   u_channel.close();
   v_channel.close();
   yuv420_channel.close();
   rgb420_channel.close();
}

int main(int argc,char**argv)
{
   std::fstream raw(argv[1],std::ios::binary | std::ios::in);

   unsigned char image[RAW_WIDTH][RAW_HEIGHT][RAW_CHANNEL];

   unsigned char data;
   for(int i = 0;i < RAW_WIDTH;i++){
      for(int j = 0;j < RAW_HEIGHT;j++){
	 for(int k = 0;k < RAW_CHANNEL;k++){
	    raw.read((char*) &data, sizeof(char));
	    image[i][j][k] = data; 
	 } 
      }
   }

   RGBChannelSplit(image,argv[2],argv[3],argv[4]);
   Y(image,argv[5]);
   I(image,argv[6]);
   Q(image,argv[7]);
   U(image,argv[8]);
   V(image,argv[9]);
   YUV420ToRGBandPSNR(image,argv[5],argv[8],argv[9]);

   raw.close();
}
