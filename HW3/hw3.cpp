/*
 *	2015/11/05
 *	Tony Guo
 *	R04458006
 *
 *	HW3 - DCT IDCT Quantization InverseQuantization RLD RLC HuffmandVLC HuffmanVLD ArithmeticVLC ArithmeticVLD
 *
 */

#include "compression.h"

enum{
   Q1 = 1,
   Q2,
   Q3,
   Q4
};

int main(int argc,char**argv)
{
   /*Q1*/
   switch(std::atoi(argv[5])){
      case Q1:
	 {
	    Compression compression(argv[1],std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
	    compression.OneD_Block_DCT();
	    compression.OneD_Block_IDCT();
	    std::cout << "Block 8x8 DCT IDCT" << std::endl;
	    std::cout << "PSNR：" << compression.PSNRComputing() << " dB" << std::endl;
	    break;
	 }

   /*Q2*/
      case Q2:
	 {
	    Compression compression(argv[1],std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
	    compression.OneD_Block_DCT();
	    compression.Quantization();
	    compression.RLC();
	    compression.RLD();
	    compression.IQuantization();
	    compression.OneD_Block_IDCT();
	    std::cout << "Block 8x8 DCT IDCT RLC RLD" << std::endl;
	    std::cout << "PSNR：" << compression.PSNRComputing() << " dB" << std::endl;
	    break;
	 }
   
   /*Q3*/
      case Q3:
	 {
	    Compression compression(argv[1],std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
	    compression.OneD_Block_DCT();
	    compression.Quantization();
	    compression.RLC();
	    compression.HuffmanVLC();
	    compression.HuffmanVLD();
	    compression.RLD();
	    compression.IQuantization();
	    compression.OneD_Block_IDCT();
	    std::cout << "Block 8x8 DCT IDCT RLC RLD HuffmanVLC HuffmanVLD" << std::endl;
	    std::cout << "CR：" << (std::stoi(argv[2]) * std::stoi(argv[3])) / compression.GetHuffmanByte() << std::endl;
	    std::cout << "PSNR：" << compression.PSNRComputing() << " dB" << std::endl;
	    break;
	 }

   /*Q4*/
      case Q4:
	 {
	    Compression compression(argv[1],std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
	    compression.OneD_Block_DCT();
	    compression.Quantization();
	    compression.RLC();
	    compression.AriVLC();
	    compression.AriVLD();
	    compression.RLD();
	    compression.IQuantization();
	    compression.OneD_Block_IDCT();
	    std::cout << "Block 8x8 DCT IDCT RLC RLD ArithmeticVLC ArithmeticVLD" << std::endl;
	    std::cout << "CR：" << (std::stoi(argv[2]) * std::stoi(argv[3])) / compression.GetAriByte() << std::endl;
	    std::cout << "PSNR：" << compression.PSNRComputing() << " dB" << std::endl;
	    break;
	 }
   }
}
