/*
 *	2015/11/22
 *	Tony Guo
 *	R04458006
 *
 *	HW4
 *
 */

#include "compression.h"

namespace{
   /*
      const double Q[8][8] = 
      {{16,  11,  10,  16,  24,  40,  51,  61}, 
	 {12,  12,  14,  19,  26,  58,  60,  55}, 
	 {14,  13,  16,  24,  40,  57,  69,  56}, 
	 {14,  17,  22,  29,  51,  87,  80,  62}, 
	 {18,  22,  37,  56,  68, 109, 103,  77}, 
	 {24,  35,  55,  64,  81, 104, 113,  92}, 
	 {49,  64,  78,  87, 103, 121, 120, 101}, 
	 {72,  92,  95,  98, 112, 100, 103,  99}};
	 */

      const int Q_INTRA[8][8] = 
      {{8,  16,  19,  22,  26,  27,  29,  34}, 
	{16,  16,  22,  24,  27,  29,  34,  37}, 
	{19,  22,  26,  27,  29,  34,  34,  38}, 
	{22,  22,  26,  27,  29,  34,  37,  40}, 
	{22,  26,  27,  29,  32,  35,  40,  48}, 
	{26,  27,  29,  32,  35,  40,  48,  58}, 
	{26,  27,  29,  34,  38,  46,  56,  69}, 
	{27,  29,  35,  38,  46,  56,  69,  83}};
  
      const int Q_INTER[8][8] = 
      {{16,  16,  16,  16,  16,  16,  16,  16}, 
	 {16,  16,  16,  16,  16,  16,  16,  16}, 
	 {16,  16,  16,  16,  16,  16,  16,  16}, 
	 {16,  16,  16,  16,  16,  16,  16,  16},
	 {16,  16,  16,  16,  16,  16,  16,  16}, 
	 {16,  16,  16,  16,  16,  16,  16,  16},
	 {16,  16,  16,  16,  16,  16,  16,  16}, 
	 {16,  16,  16,  16,  16,  16,  16,  16}};
}		  	  

void Compression::LoadRawImage(std::vector<std::vector<unsigned char> > &origin_image)
{
   if(origin_image.size()){
      origin_image.clear();
   }

   std::fstream raw(path,std::ios::binary | std::ios::in);

   for(int x = 0;x < raw_width;++x){
      std::vector<unsigned char> temp;
      for(int y = 0;y < raw_height;++y){
	 unsigned char data;
	 raw.read((char*) &data, sizeof(char));
	 temp.push_back(data);
      }
      origin_image.push_back(temp);
   }

   raw.close();
}

void Compression::LoadRawVideo(std::vector<std::vector<std::vector<unsigned char> > > &video)
{
   if(video.size()){
      video.clear();
   }

   std::fstream raw(path,std::ios::binary | std::ios::in);
   std::fstream raw_end(path,std::ios::binary | std::ios::in | std::ios::ate);

   total_frame = raw_end.tellg() / (raw_width * raw_height);
   for(int f = 0;f < raw_end.tellg() / (raw_width * raw_height);++f){
      std::vector<std::vector<unsigned char> > image;
      for(int x = 0;x < raw_width;++x){
	 std::vector<unsigned char> temp;
	 for(int y = 0;y < raw_height;++y){
	    unsigned char data;
	    raw.read((char*) &data, sizeof(char));
	    temp.push_back(data);
	 }
	 image.push_back(temp);
      }
      video.push_back(image);
   }

   raw.close();
   raw_end.close();
}

void Compression::OneD_Block_DCT(std::vector<std::vector<unsigned char> > &origin_image,std::vector<std::vector<double> > &dct)
{
   if(dct.size()){
      dct.clear();
   }

   for(int u = 0;u < raw_width;++u){
      std::vector<double> temp;
      for(int v = 0;v < raw_height;++v){
	 temp.push_back(0.0f);
      }
      dct.push_back(temp);
   }

   double dct_temp[raw_width][raw_height];
   for(int i = 0;i < raw_width / block_size;++i){
      for(int j = 0;j < raw_height / block_size;++j){
	 for(int u = 0;u < block_size;++u){
	    for(int v = 0;v < block_size;++v){
	       double f_uv = 0;
	       for(int y = 0;y < block_size;++y){
		  f_uv += (double)origin_image[i * block_size + u][j * block_size + y] * std::cos((M_PI * (2 * y + 1) * v) / (2 * block_size));
	       }
	       dct_temp[i * block_size + u][j * block_size + v] = f_uv * C(v) * std::sqrt(2.0f / block_size);
	    }
	 }
      }
   }

   for(int i = 0;i < raw_width / block_size;++i){
      for(int j = 0;j < raw_height / block_size;++j){
	 for(int v = 0;v < block_size;++v){
	    for(int u = 0;u < block_size;++u){
	       double f_uv = 0;
	       for(int x = 0;x < block_size;++x){
		  f_uv += (double)dct_temp[i * block_size + x][j * block_size + v] * std::cos((M_PI * (2 * x + 1) * u) / (2 * block_size));
	       }
	       dct[i * block_size + u][j * block_size + v] = f_uv * C(u) * std::sqrt(2.0f / block_size);
	    }
	 }
      }
   }

   /*
   std::fstream raw(path + std::string("_dct.raw"),std::ios::binary | std::ios::out);
   for(int i = 0;i < dct.size();++i){
      for(int j = 0;j < dct[i].size();++j){
	 raw << (unsigned char)limit(dct[i][j] + 128.0f);
      }
   }
   */
}

void Compression::OneD_Block_IDCT(std::vector<std::vector<double> > &dct,std::vector<std::vector<unsigned char> > &idct)
{
   if(idct.size()){
      idct.clear();
   }

   for(int u = 0;u < raw_width;++u){
      std::vector<unsigned char> temp;
      for(int v = 0;v < raw_height;++v){
	 temp.push_back(0.0f);
      }
      idct.push_back(temp);
   }

   double idct_temp[raw_width][raw_height];
   for(int i = 0;i < raw_width / block_size;++i){
      for(int j = 0;j < raw_height / block_size;++j){
	 for(int x = 0;x < block_size;++x){
	    for(int y = 0;y < block_size;++y){
	       double f_xy = 0;
	       for(int v = 0;v < block_size;++v){
		  f_xy += (double)dct[i * block_size + x][j * block_size + v] * C(v) * std::cos((M_PI * (2 * y + 1) * v) / (2 * block_size));
	       }
	       idct_temp[i * block_size + x][j * block_size + y] = f_xy * std::sqrt(2.0f / block_size);
	    }
	 }
      }
   }

   for(int i = 0;i < raw_width / block_size;++i){
      for(int j = 0;j < raw_height / block_size;++j){
	 for(int y = 0;y < block_size;++y){
	    for(int x = 0;x < block_size;++x){
	       double f_xy = 0;
	       for(int u = 0;u < block_size;++u){
		  f_xy += (double)idct_temp[i * block_size + u][j * block_size + y] * C(u) * std::cos((M_PI * (2 * x + 1) * u) / (2 * block_size));
	       }
	       f_xy *= std::sqrt(2.0f / block_size);
	       if(f_xy > 255){idct[i * block_size + x][j * block_size + y] = 255;}
	       else if(f_xy < 0){idct[i * block_size + x][j * block_size + y] = 0;}
	       else{idct[i * block_size + x][j * block_size + y] = (unsigned char)std::round(f_xy);}
	    }
	 }
      }
   }

   
   /*
   std::fstream raw(path + std::string("_idct.raw"),std::ios::binary | std::ios::out | std::ios::app);
   for(int i = 0;i < dct.size();++i){
      for(int j = 0;j < dct[i].size();++j){
	 raw << idct[i][j];
      }
   }
   */
   
}

void Compression::Quantization(std::vector<std::vector<double> > &dct,const int mode)
{
   const int (*Q)[8] = (mode == INTRA)?&(Q_INTRA[0]):&(Q_INTER[0]);
   for(int x = 0;x < raw_width;x += block_size){
      for(int y = 0;y < raw_height;y += block_size){
	 for(int u = 0;u < block_size;++u){
	    for(int v = 0;v < block_size;++v){
	       dct[x + u][y + v] = (int)std::round(dct[x + u][y + v] / Q[u][v]);
	    }
	 }
      }
   }
}

void Compression::IQuantization(std::vector<std::vector<double> > &dct,const int mode)
{
   const int (*Q)[8] = (mode == INTRA)?&(Q_INTRA[0]):&(Q_INTER[0]);
   for(int x = 0;x < raw_width;x += block_size){
      for(int y = 0;y < raw_height;y += block_size){
	 for(int u = 0;u < block_size;++u){
	    for(int v = 0;v < block_size;++v){
	       dct[x + u][y + v] *= Q[u][v];
	    }
	 }
      }
   }
}

int idx;

void Compression::ME(std::vector<std::vector<unsigned char> > &origin_image,std::vector<std::vector<unsigned char> > &frame_buffer,int d,unsigned int mode)
{
   int SAD, x, y;
   if(mode == FULL){
      for(int m = 0; m < raw_width; m += block_size){
	 for(int n = 0; n < raw_height; n += block_size){
	    std::map<int, std::pair<int,int> > diff;
	    for(int i = -d; i < d; ++i){
	       for(int j = -d; j < d; ++j){
		  if((m + i < 0) || (m + i + block_size > raw_width - 1) || (n + j < 0) || (n + j + block_size > raw_height - 1)){
		  }
		  else{
		     SAD = 0;
		     for(int k = 0; k < block_size; ++k){
			for(int l = 0; l < block_size; ++l){
			   SAD += std::abs((int)origin_image[m + k][n + l] - (int)frame_buffer[m + i + k][n + j + l]);
			}
		     }
		     diff[SAD] = std::make_pair(m + i, n + j);
		  }
	       }
	    }
	    MV.push_back(diff.begin()->second);
	 }
      }
   }
   else if(mode == FAST){
      for(int m = 0; m < raw_width; m += block_size){
	 for(int n = 0; n < raw_height; n += block_size){
	    x = y = -4;
	    for(int step = 4; step > 0; step/=2){
	       std::map<int, std::pair<int,int> > diff;
	       for(int i = x; i < x + 2 * step + 1; i += step){
		  for(int j = y; j < y + 2 * step + 1; j += step){
		  if((m + i < 0) || (m + i + block_size > raw_width - 1) || (n + j < 0) || (n + j + block_size > raw_height - 1)){
		  }
		  else{
			SAD = 0;
			for(int k = 0; k < block_size; ++k){
			   for(int l = 0; l < block_size; ++l){
			      SAD += std::abs((int)origin_image[m + k][n + l] - (int)frame_buffer[m + i + k][n + j + l]);
			   }
			}
			diff[SAD] = std::make_pair(i,j);
		     }
		  }
	       }	
	       x = diff.begin()->second.first - step / 2;
	       y = diff.begin()->second.second - step / 2;
	       if(step/2 == 0) MV.push_back(std::make_pair(diff.begin()->second.first + m, diff.begin()->second.second + n));	
	    }
	 }
      }
   }
}

void Compression::MC(std::vector<std::vector<unsigned char> > &frame_buffer,std::vector<std::vector<unsigned char> > &predict_frame)
{
   idx = 0;
   for(int m = 0; m < raw_width; m += block_size){
      for(int n = 0; n < raw_height; n += block_size){
	 for(int k = 0; k < block_size; ++k){
	    for(int l = 0; l < block_size; ++l){
	       predict_frame[m + k][n + l] = frame_buffer[MV[idx].first + k][MV[idx].second + l];
	    }
	 }
	 ++idx;
      }
   }
   MV.clear();
}

void Compression::NextFrame(std::vector<std::vector<std::vector<unsigned char> > > &video,std::vector<std::vector<unsigned char> > &origin_image,unsigned int frame_num)
{
   if(origin_image.size()){
      origin_image.clear();
   }

   for(int x = 0;x < raw_width;++x){
      std::vector<unsigned char> temp;
      for(int y = 0;y < raw_height;++y){
	 temp.push_back(video[frame_num][x][y]);
      }
      origin_image.push_back(temp);
   }
}
