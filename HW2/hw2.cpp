/*
 *	2015/10/18
 *	Tony Guo
 *	R04458006
 *
 *	HW2 - DCT IDCT Quantization
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

const int RAW_WIDTH = 64;
const int RAW_HEIGHT = 64;
const int N = 64;

inline double C(double value)
{
  return (value == 0)?1/std::sqrt(2):1.0f;
}

inline double limit(double value)
{
   return (value < 0)?0:(value > 255)?255:value;
}

inline int Truncate_Shift(double value,int bits)
{
   value = std::floor(value);
   int shift = 0;
   while(value >= (2 << ((bits-2))) || value < -(2 << ((bits-2)))){
      value /= 2.0f;
      ++shift;
   }

   return shift;
}

double PSNRComputing(unsigned char (&origi_image)[RAW_WIDTH][RAW_HEIGHT],unsigned char (&trans_image)[RAW_WIDTH][RAW_HEIGHT])
{
   double mse = 0,psnr = 0;

   for(int i = 0;i < RAW_WIDTH;i++){
      for(int j = 0;j < RAW_HEIGHT;j++){
	 mse += (origi_image[i][j] - trans_image[i][j]) * (origi_image[i][j] - trans_image[i][j]);
      }
   }

   if(mse != 0){
      mse /= RAW_WIDTH * RAW_HEIGHT;
      psnr = 10 * std::log10((255 * 255) / mse);
   }
   else{
      psnr = 999999;
   }

   //std::cout << "PSNR : " << psnr << " dB" << std::endl;
   return psnr;
}

void TwoD_Block_DCT(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT],double (&dct)[N][N],int block)
{
   for(int i = 0;i < RAW_WIDTH / block;++i){
      for(int j = 0;j < RAW_HEIGHT / block;++j){
	 for(int u = 0;u < block;++u){
	    for(int v = 0;v < block;++v){
	       double f_uv = 0;
	       for(int x = 0;x < block;++x){
		  for(int y = 0;y < block;++y){
		     f_uv += (double)image[i * block + x][j * block + y] * std::cos((M_PI * (2 * x + 1) * u) / (2 * block)) * 
			std::cos((M_PI * (2 * y + 1) * v) / (2 * block));
		  }
	       }
	       dct[i * block + u][j * block + v] = f_uv * C(u) * C(v) * (2.0f / (double)block);
	    }
	 }
      }
   }
}

void TwoD_Block_IDCT(double (&dct)[N][N],unsigned char (&idct)[RAW_WIDTH][RAW_HEIGHT],int block)
{
   for(int i = 0;i < RAW_WIDTH / block;++i){
      for(int j = 0;j < RAW_HEIGHT / block;++j){
	 for(int x = 0;x < block;++x){
	    for(int y = 0;y < block;++y){
	       double f_xy = 0;
	       for(int u = 0;u < block;++u){
		  for(int v = 0;v < block;++v){
		     f_xy += dct[i * block + u][j * block + v] * C(u) * C(v) * std::cos((M_PI * (2 * x + 1) * u) / (2 * block)) * 
			std::cos((M_PI * (2 * y + 1) * v) / (2 * block));
		  }
	       }
	       idct[i * block + x][j * block + y] = f_xy * (2.0f / (double)block);
	    }
	 }
      }
   }
}

void OneD_Block_DCT(unsigned char (&image)[RAW_WIDTH][RAW_HEIGHT],double (&dct)[N][N],int block)
{
   double dct_temp[N][N];
   for(int i = 0;i < RAW_WIDTH / block;++i){
      for(int j = 0;j < RAW_HEIGHT / block;++j){
	 for(int u = 0;u < block;++u){
	    for(int v = 0;v < block;++v){
	       double f_uv = 0;
	       for(int y = 0;y < block;++y){
		  f_uv += (double)image[i * block + u][j * block + y] * std::cos((M_PI * (2 * y + 1) * v) / (2 * block));
	       }
	       dct_temp[i * block + u][j * block + v] = f_uv * C(v) * std::sqrt(2.0f / block);
	    }
	 }
      }
   }

   for(int i = 0;i < RAW_WIDTH / block;++i){
      for(int j = 0;j < RAW_HEIGHT / block;++j){
	 for(int v = 0;v < block;++v){
	    for(int u = 0;u < block;++u){
	       double f_uv = 0;
	       for(int x = 0;x < block;++x){
		  f_uv += (double)dct_temp[i * block + x][j * block + v] * std::cos((M_PI * (2 * x + 1) * u) / (2 * block));
	       }
	       dct[i * block + u][j * block + v] = f_uv * C(u) * std::sqrt(2.0f / block);
	    }
	 }
      }
   }
}

void OneD_Block_IDCT(double (&dct)[N][N],unsigned char (&idct)[RAW_WIDTH][RAW_HEIGHT],int block)
{
   double idct_temp[N][N];
   for(int i = 0;i < RAW_WIDTH / block;++i){
      for(int j = 0;j < RAW_HEIGHT / block;++j){
	 for(int x = 0;x < block;++x){
	    for(int y = 0;y < block;++y){
	       double f_xy = 0;
	       for(int v = 0;v < block;++v){
		  f_xy += (double)dct[i * block + x][j * block + v] * C(v) * std::cos((M_PI * (2 * y + 1) * v) / (2 * block));
	       }
	       idct_temp[i * block + x][j * block + y] = f_xy * std::sqrt(2.0f / block);
	    }
	 }
      }
   }

   for(int i = 0;i < RAW_WIDTH / block;++i){
      for(int j = 0;j < RAW_HEIGHT / block;++j){
	 for(int y = 0;y < block;++y){
	    for(int x = 0;x < block;++x){
	       double f_xy = 0;
	       for(int u = 0;u < block;++u){
		  f_xy += (double)idct_temp[i * block + u][j * block + y] * C(u) * std::cos((M_PI * (2 * x + 1) * u) / (2 * block));
	       }
	       idct[i * block + x][j * block + y] = (unsigned char)std::round(f_xy * std::sqrt(2.0f / block));
	    }
	 }
      }
   }
}

void Quantization(double (&dct)[N][N],unsigned char (&idct)[RAW_WIDTH][RAW_HEIGHT],bool dead_zone_flag,int dead_zone_min,int dead_zone_max)
{
   int dc_shift = Truncate_Shift(dct[0][0],16);

   double max = 0;
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 if(u == 0 && v == 0){}
	 else{
	    if(std::fabs(dct[u][v]) > max){max = std::fabs(dct[u][v]);}
	 }
      }
   }
   int ac_shift = Truncate_Shift(max,8);

   int truncate_bits = 0;
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 if(u == 0 && v == 0){dct[u][v] = std::floor(dct[u][v]) / (1 << dc_shift);}
	 else{
	    dct[u][v] = std::floor(dct[u][v]) / (1 << ac_shift);
	    if(dead_zone_flag){
	       if(dct[u][v] <= dead_zone_max && dct[u][v] >= dead_zone_min){dct[u][v] = 0;truncate_bits += 8;}
	    }
	 }
      }
   }

   std::cout << "Dead Zone：" << dead_zone_min << " ~ " << dead_zone_max << std::endl;
   std::cout << "DC-AC Truncate：" << dc_shift << " , " << ac_shift << std::endl;
   std::cout << "Truncate Bits / Total Bits：" << truncate_bits << " / " << 64 * 64 * 8 - 1 << " = " << (double)truncate_bits / (64 * 64 * 8 - 1) << std::endl;

   double idct_temp[N][N];
   for(int x = 0;x < RAW_WIDTH;++x){
      for(int y = 0;y < RAW_HEIGHT;++y){
	 double f_xy = 0;
	 for(int v = 0;v < N;++v){
	    if(x == 0 && v == 0){f_xy += (double)((int)dct[x][v] << dc_shift) * C(v) * std::cos((M_PI * (2 * y + 1) * v) / (2 * N));}
	    else{f_xy += (double)((int)dct[x][v] << ac_shift) * C(v) * std::cos((M_PI * (2 * y + 1) * v) / (2 * N));}
	 }
	 idct_temp[x][y] = f_xy * std::sqrt(2.0f / N);
      }
   }
   for(int y = 0;y < RAW_HEIGHT;++y){
      for(int x = 0;x < RAW_WIDTH;++x){
	 double f_xy = 0;
	 for(int u = 0;u < N;++u){
	    if(u == 0 && y == 0){f_xy += (double)idct_temp[u][y] * C(u) * std::cos((M_PI * (2 * x + 1) * u) / (2 * N));}
	    else{f_xy += (double)idct_temp[u][y] * C(u) * std::cos((M_PI * (2 * x + 1) * u) / (2 * N));}
	 }
	 idct[x][y] = (unsigned char)limit(std::round(f_xy * std::sqrt(2.0f / N)));
      }
   }

}

int main(int argc,char**argv)
{
   std::fstream raw(argv[1],std::ios::binary | std::ios::in);

   unsigned char image[RAW_WIDTH][RAW_HEIGHT];

   unsigned char data;
   for(int i = 0;i < RAW_WIDTH;i++){
      for(int j = 0;j < RAW_HEIGHT;j++){
	 raw.read((char*) &data, sizeof(char));
	 image[i][j] = data; 
      }
   }

   int block = 64;
   
   double dct[N][N];
   clock_t dct_time = std::clock();
   TwoD_Block_DCT(image,dct,block);
   std::cout << "2D-DCT：" << (double)(clock() - dct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream two_d_dct_file(std::string(argv[1]) + std::string("_2d_dct.raw"),std::ios::binary | std::ios::out);
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 two_d_dct_file << (unsigned char)limit((dct[u][v] + 128.0f));
      }
   }

   unsigned char idct[RAW_WIDTH][RAW_HEIGHT];
   clock_t idct_time = std::clock();
   TwoD_Block_IDCT(dct,idct,block);
   std::cout << "2D-IDCT：" << (double)(clock() - idct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream two_d_idct_file(std::string(argv[1]) + std::string("_2d_idct.raw"),std::ios::binary | std::ios::out);
   for(int x = 0;x < N;++x){
      for(int y = 0;y < N;++y){
	 two_d_idct_file << idct[x][y];
      }
   }

   std::cout << "2D-IDCT PSNR：" << PSNRComputing(image,idct) << std::endl;

   dct_time = std::clock();
   OneD_Block_DCT(image,dct,block);
   std::cout << "1D-DCT：" << (double)(clock() - dct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream one_d_dct_file(std::string(argv[1]) + std::string("_1d_dct.raw"),std::ios::binary | std::ios::out);
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 one_d_dct_file << (unsigned char)limit((dct[u][v] + 128.0f));
      }
   }

   idct_time = std::clock();
   OneD_Block_IDCT(dct,idct,block);
   std::cout << "1D-IDCT：" << (double)(clock() - idct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream one_d_idct_file(std::string(argv[1]) + std::string("_1d_idct.raw"),std::ios::binary | std::ios::out);
   for(int x = 0;x < N;++x){
      for(int y = 0;y < N;++y){
	 one_d_idct_file << idct[x][y];
      }
   }

   std::cout << "1D-IDCT PSNR：" << PSNRComputing(image,idct) << std::endl;

   double quantization_0_0[N][N];
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 quantization_0_0[u][v] = dct[u][v];
      }
   }

   Quantization(quantization_0_0,idct,false,0,0);
   std::fstream quantization_0_0_idct_file(std::string(argv[1]) + std::string("_quantization_0_0_idct.raw"),std::ios::binary | std::ios::out);
   for(int x = 0;x < N;++x){
      for(int y = 0;y < N;++y){
	 quantization_0_0_idct_file << idct[x][y];
      }
   }

   std::cout << "Quantization PSNR：" << PSNRComputing(image,idct) << std::endl;

   double quantization__7_7[N][N];
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 quantization__7_7[u][v] = dct[u][v];
      }
   }


   Quantization(quantization__7_7,idct,true,-7,7);
   std::fstream quantization__7_7_idct_file(std::string(argv[1]) + std::string("_quantization__7_7_idct.raw"),std::ios::binary | std::ios::out);
   for(int x = 0;x < N;++x){
      for(int y = 0;y < N;++y){
	 quantization__7_7_idct_file << idct[x][y];
      }
   }

   std::cout << "Quantization Truncate -7~7 PSNR：" << PSNRComputing(image,idct) << std::endl;

   block = 8;

   dct_time = std::clock();
   TwoD_Block_DCT(image,dct,block);
   std::cout << "2D-Block-DCT：" << (double)(clock() - dct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream two_d_block_dct_file(std::string(argv[1]) + std::string("_2d_block_dct.raw"),std::ios::binary | std::ios::out);
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 two_d_block_dct_file << (unsigned char)limit((dct[u][v] + 128.0f));
      }
   }

   idct_time = std::clock();
   TwoD_Block_IDCT(dct,idct,block);
   std::cout << "2D-Block-IDCT：" << (double)(clock() - idct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream two_d_block_idct_file(std::string(argv[1]) + std::string("_2d_block_idct.raw"),std::ios::binary | std::ios::out);
   for(int x = 0;x < N;++x){
      for(int y = 0;y < N;++y){
	 two_d_block_idct_file << idct[x][y];
      }
   }

   std::cout << "2D-Block-IDCT PSNR：" << PSNRComputing(image,idct) << std::endl;

   dct_time = std::clock();
   OneD_Block_DCT(image,dct,block);
   std::cout << "1D-Block-DCT：" << (double)(clock() - dct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream one_d_block_dct_file(std::string(argv[1]) + std::string("_1d_block_dct.raw"),std::ios::binary | std::ios::out);
   for(int u = 0;u < N;++u){
      for(int v = 0;v < N;++v){
	 one_d_block_dct_file << (unsigned char)limit((dct[u][v] + 128.0f));
      }
   }

   idct_time = std::clock();
   OneD_Block_IDCT(dct,idct,block);
   std::cout << "1D-Block-IDCT：" << (double)(clock() - idct_time) / CLOCKS_PER_SEC << "s" << std::endl;

   std::fstream one_d_block_idct_file(std::string(argv[1]) + std::string("_1d_block_idct.raw"),std::ios::binary | std::ios::out);
   for(int x = 0;x < N;++x){
      for(int y = 0;y < N;++y){
	 one_d_block_idct_file << idct[x][y];
      }
   }

   std::cout << "1D-Block-IDCT PSNR：" << PSNRComputing(image,idct) << std::endl;

   raw.close();
   two_d_dct_file.close();
   two_d_idct_file.close();
   one_d_dct_file.close();
   one_d_idct_file.close();
   quantization_0_0_idct_file.close();
   quantization__7_7_idct_file.close();
   two_d_block_dct_file.close();
   two_d_block_idct_file.close();
   one_d_block_dct_file.close();
   one_d_block_idct_file.close();
}
