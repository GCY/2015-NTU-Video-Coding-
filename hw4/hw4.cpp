/*
 *	2015/11/22
 *	Tony Guo
 *	R04458006
 *
 *	HW4
 *
 */

#include "compression.h"

const int GOP = 16;

void WriteFrame2File(std::fstream &file,std::vector<std::vector<unsigned char> > frame)
{
   for(int i = 0;i < frame.size();++i){
      for(int j = 0;j < frame[i].size();++j){
	 file << frame[i][j];
      }
   }
}

int main(int argc,char**argv)
{
   std::string path = argv[1];
   unsigned int width = std::stoi(argv[2]);
   unsigned int height = std::stoi(argv[3]);
   unsigned int block_size = std::stoi(argv[4]);
   int d = std::stoi(argv[5]);
   unsigned int mode = std::stoi(argv[6]);

   std::stringstream ds;
   ds << d;
   std::fstream diff_file(path + std::string((mode == FULL)?".full":".fast") + std::string(".W") + ds.str() + std::string(".diff"),std::ios::out | std::ios::binary | std::ios::app);
   std::fstream res_file(path + std::string((mode == FULL)?".full":".fast") + std::string(".W") + ds.str() + std::string(".res"),std::ios::out | std::ios::binary | std::ios::app);
   std::fstream rec_file(path + std::string((mode == FULL)?".full":".fast") + std::string(".W") + ds.str() + std::string(".rec"),std::ios::out | std::ios::binary | std::ios::app);

   std::fstream time_file(path + std::string((mode == FULL)?".full":".fast") + std::string(".W") + ds.str() + std::string(".time.txt"),std::ios::out);
   std::fstream psnr_file(path + std::string((mode == FULL)?".full":".fast") + std::string(".W") + ds.str() + std::string(".psnr.txt"),std::ios::out);


   std::vector<std::vector<std::vector<unsigned char> > > video;

   std::vector<std::vector<unsigned char> > frame_buffer;
   std::vector<std::vector<unsigned char> > predict_frame;
   std::vector<std::vector<unsigned char> > diff_frame;
   std::vector<std::vector<unsigned char> > res_frame;
   for(int i = 0;i < width;++i){
      std::vector<unsigned char> temp;
      for(int j = 0;j < height;++j){
	 temp.push_back(0);
      }
      frame_buffer.push_back(temp);
      predict_frame.push_back(temp);
      diff_frame.push_back(temp);
      res_frame.push_back(temp);
   }

   Compression compression(path,width,height,block_size);
   compression.LoadRawVideo(video);

   for(int i = 0;i < compression.TotalFrame();++i){
      std::cout << "frame" << i << std::endl;
      clock_t time = std::clock();
      std::vector<std::vector<unsigned char> > origin_image;
      std::vector<std::vector<double> > dct;
      compression.NextFrame(video,origin_image,i);
      if((i % GOP) == INTRA){
	 WriteFrame2File(diff_file,frame_buffer);
	 compression.OneD_Block_DCT(origin_image,dct);
	 compression.Quantization(dct,INTRA);
	 compression.IQuantization(dct,INTRA);
	 compression.OneD_Block_IDCT(dct,frame_buffer);
	 WriteFrame2File(res_file,frame_buffer);
	 double PSNR = compression.PSNRComputing(origin_image,frame_buffer);
	 std::cout << "PSNR：" << compression.PSNRComputing(origin_image,frame_buffer) << " dB" << std::endl;
	 psnr_file << PSNR << std::endl;
      }
      else{
	 compression.ME(origin_image,frame_buffer,d,mode);
	 compression.MC(frame_buffer,predict_frame);
	 WriteFrame2File(rec_file,predict_frame);
	 compression.SubImage(origin_image,predict_frame,diff_frame);
	 WriteFrame2File(diff_file,diff_frame);
	 compression.OneD_Block_DCT(diff_frame,dct);
	 compression.Quantization(dct,INTRA);
	 compression.IQuantization(dct,INTRA);
	 compression.OneD_Block_IDCT(dct,res_frame);
	 WriteFrame2File(res_file,res_frame);

	 compression.AddImage(res_frame,predict_frame,frame_buffer);
	 double PSNR = compression.PSNRComputing(origin_image,frame_buffer);
	 std::cout << "PSNR：" << PSNR << " dB" << std::endl;
	 psnr_file << PSNR << std::endl;
      }
      double diff_time = (double)(clock() - time) / CLOCKS_PER_SEC;
      std::cout << "time：" << diff_time << "s" << std::endl;
      time_file << diff_time << std::endl;
   }

   diff_file.close();
   res_file.close();
   rec_file.close();

   time_file.close();
   psnr_file.close();

}
