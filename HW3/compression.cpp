/*
 *	2015/11/05
 *	Tony Guo
 *	R04458006
 *
 *	HW3 - DCT IDCT Quantization InverseQuantization RLD RLC HuffmandVLC HuffmanVLD ArithmeticVLC ArithmeticVLD
 *
 */

#include "compression.h"

namespace{
      const double Q[8][8] = 
      {{16,  11,  10,  16,  24,  40,  51,  61}, 
	 {12,  12,  14,  19,  26,  58,  60,  55}, 
	 {14,  13,  16,  24,  40,  57,  69,  56}, 
	 {14,  17,  22,  29,  51,  87,  80,  62}, 
	 {18,  22,  37,  56,  68, 109, 103,  77}, 
	 {24,  35,  55,  64,  81, 104, 113,  92}, 
	 {49,  64,  78,  87, 103, 121, 120, 101}, 
	 {72,  92,  95,  98, 112, 100, 103,  99}};
		  	  
      const int zigzag_x[8][8] = 		  	  
      {{0, 0, 1, 2, 1, 0, 0, 1}, 
	 {2, 3, 4, 3, 2, 1, 0, 0}, 
	 {1, 2, 3, 4, 5, 6, 5, 4}, 
	 {3, 2, 1, 0, 0, 1, 2, 3}, 
	 {4, 5, 6, 7, 7, 6, 5, 4}, 
	 {3, 2, 1, 2, 3, 4, 5, 6}, 
	 {7, 7, 6, 5, 4, 3, 4, 5}, 
	 {6, 7, 7, 6, 5, 6, 7, 7}};
  
      const int zigzag_y[8][8] = 
      {{0, 1, 0, 0, 1, 2, 3, 2}, 
	 {1, 0, 0, 1, 2, 3, 4, 5}, 
	 {4, 3, 2, 1, 0, 0, 1, 2}, 
	 {3, 4, 5, 6, 7, 6, 5, 4}, 
	 {3, 2, 1, 0, 1, 2, 3, 4}, 
	 {5, 6, 7, 7, 6, 5, 4, 3}, 
	 {2, 3, 4, 5, 6, 7, 7, 6}, 
	 {5, 4, 5, 6, 7, 7, 6, 7}};
}

void Compression::LoadRawImage()
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

void Compression::OneD_Block_DCT()
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

   std::fstream raw(path + std::string("_dct.raw"),std::ios::binary | std::ios::out);
   for(int i = 0;i < dct.size();++i){
      for(int j = 0;j < dct[i].size();++j){
	 raw << (unsigned char)limit(dct[i][j] + 128.0f);
      }
   }
}

void Compression::OneD_Block_IDCT()
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

   std::fstream raw(path + std::string("_idct.raw"),std::ios::binary | std::ios::out);
   for(int i = 0;i < dct.size();++i){
      for(int j = 0;j < dct[i].size();++j){
	 raw << idct[i][j];
      }
   }
}

void Compression::Quantization()
{
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

void Compression::IQuantization()
{
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

void Compression::RLC()
{
   std::fstream rlc(path + std::string(".rlc"),std::ios::out);

   for(int x = 0;x < raw_width;x += block_size){
      for(int y = 0;y < raw_height;y += block_size){
	 unsigned int zero_count = 0;
	 for(int u = 0;u < block_size;++u){
	    for(int v = 0;v < block_size;++v){
	       if(dct[x + zigzag_x[u][v]][y + zigzag_y[u][v]]){
		  rlc << "(" << zero_count << "," << dct[x + zigzag_x[u][v]][y + zigzag_y[u][v]] << ")";
		  zero_count = 0;
	       }
	       else{
		  ++zero_count;
	       }
	    }
	 }
	 if(zero_count){
	    rlc << "EOB";
	    rlc << std::endl;
	 }
	 else{
	    rlc << std::endl;
	 }
      }
   }

   rlc.close();
}

void Compression::RLD()
{
   std::fstream rld(path + std::string(".rlc"),std::ios::in);

   std::vector<std::string> rlc;
   std::vector<std::vector<int> > rld_table;
   std::vector<std::vector<int> > number_table;

   std::string str;
   while(rld >> str){
      rlc.push_back(str);
   }

   for(size_t i = 0;i < rlc.size();++i){
      std::vector<int> temp;
      std::string number;
      for(size_t j = 0;j < rlc[i].size();++j){
	 if((rlc[i][j] >= '0' && rlc[i][j] <= '9') || rlc[i][j] == '-'){
	    std::stringstream ss;
	    ss << rlc[i][j];
	    number += ss.str();
	 }
	 else{
	    if(number.size()){
	       temp.push_back(std::stoi(number));
	       number.clear();
	    }
	 }
      }
      number_table.push_back(temp);
   }

   for(size_t i = 0;i < number_table.size();++i){
      std::vector<int> temp;
      for(size_t j = 0;j < number_table[i].size();++j){
	 if(!(j % 2)){
	    for(int k = number_table[i][j];k > 0;--k){
	       temp.push_back(0);
	    }
	 }
	 else{
	    if(number_table[i][j] != 0){
	       temp.push_back(number_table[i][j]);
	    }
	 }
      }
      for(size_t j = temp.size();j < (block_size * block_size);++j){
	 temp.push_back(0);
      }
      rld_table.push_back(temp);
   }

   for(int x = 0;x < raw_width;x += block_size){
      for(int y = 0;y < raw_height;y += block_size){
	 unsigned int count = 0;
	 for(size_t i = 0;i < rld_table[x + (y / block_size)].size();++i){
	    dct[x + zigzag_x[count / block_size][count % block_size]][y + zigzag_y[count / block_size][count % block_size]] = rld_table[x + (y / block_size)][i];
	    ++count;
	 }
      }
   }


   rld.close();

}

void Compression::CreateHuffmanTree(std::map<char,int> &prob_table,std::vector<HuffmanNode*> &huffman_array)
{
   for(std::map<char,int>::iterator it = prob_table.begin();it != prob_table.end();++it){
      huffman_array.push_back(new HuffmanNode(it->second,it->first,std::string(""),NULL,NULL));
   }

   for(size_t i = 0;i < huffman_array.size() - 1;++i){
      std::sort(huffman_array.begin(),huffman_array.end(),HuffmanArraySort);
      HuffmanNode *node1 = huffman_array[i];
      HuffmanNode *node2 = huffman_array[i+1];
      HuffmanNode *new_node = new HuffmanNode(node1->prob + node2->prob,0,std::string(),node1,node2);
      huffman_array[i+1] = new_node;
   }
}

void Compression::CreateHuffmanCode(std::vector<HuffmanNode> &huffman_table,HuffmanNode *root,std::string code)
{
   root->code = code;
   if(root->data != 0){
      huffman_table.push_back(HuffmanNode(root->prob,root->data,root->code,NULL,NULL));
   }
   if(root->L){CreateHuffmanCode(huffman_table,root->L,root->code + std::string("0"));}
   if(root->R){CreateHuffmanCode(huffman_table,root->R,root->code + std::string("1"));}
}

void Compression::HuffmanVLC()
{
   std::fstream rld(path + std::string(".rlc"),std::ios::in);
   std::fstream huf(path + std::string(".huf"),std::ios::out);

   std::vector<std::string> rlc;
   std::string str;
   while(rld >> str){
      rlc.push_back(str);
   }

   std::map<char,int> prob_table;
   for(size_t i = 0;i < rlc.size();++i){
      for(size_t j = 0;j < rlc[i].size();++j){
	 if(prob_table.insert(std::make_pair(rlc[i][j],1)).second){
	 }
	 else{
	    ++prob_table[rlc[i][j]];
	 }
      }
   }

   for(std::map<char,int>::iterator it = prob_table.begin();it != prob_table.end();++it){
      huf << it->first << "/" << it->second << "/";
   }
   huf << std::endl;

   std::vector<HuffmanNode*> huffman_array;
   CreateHuffmanTree(prob_table,huffman_array);

   std::vector<HuffmanNode> huffman_table;
   CreateHuffmanCode(huffman_table,huffman_array[huffman_array.size()-1],std::string(""));

   for(size_t i = 0;i < rlc.size();++i){
      for(size_t j = 0;j < rlc[i].size();++j){
	 huf << std::find(huffman_table.begin(),huffman_table.end(),rlc[i][j])->code;
	 huffman_total_byte += std::find(huffman_table.begin(),huffman_table.end(),rlc[i][j])->code.size();
      }
      huf << std::endl;
   }

   huffman_total_byte /= 8.0f;
   
/* //binary file
   std::vector<std::vector<bool> > binary;
   for(size_t i = 0;i < rlc.size();++i){
      std::vector<bool> temp;
      for(size_t j = 0;j < rlc[i].size();++j){
	 std::string code = std::find(huffman_table.begin(),huffman_table.end(),rlc[i][j])->code;
	 for(size_t k = 0;k < code.size();++k){
	    temp.push_back(((code[k] == '0')?0:1));
	 }
      }
      binary.push_back(temp);
   }
   
   char c;
   for(size_t i = 0;i < binary.size();++i){
      for(size_t j = 0;j < binary[i].size();++j){
	 if((j % 8) == 0){
	    huf << c;
	    c = 0;
	 }
	 c += (binary[i][j] << (j % 8));
      }
   }

   //

   std::fstream huf2(path + std::string(".huf"),std::ios::in);
   std::fstream hufd(path + std::string(".ddhuf"),std::ios::out);
   huf2 >> str;
   hufd << str;
   hufd << std::endl;
   while(!huf2.eof()){
      huf2.get(c);
      for(int j = 0;j < 8;++j){
	 hufd << ((c & (1 << j))?'1':'0');
      }
      hufd << std::endl;
   }
*/
   rld.close();
   huf.close();
}

void Compression::HuffmanVLD()
{
   std::fstream huf(path + std::string(".huf"),std::ios::in);
   std::fstream rld(path + std::string(".rlc"),std::ios::out);

   std::string str;
   huf >> str;
   std::stringstream ss(str);

   std::vector<std::string> split;
   for (std::string each;std::getline(ss, each, '/');split.push_back(each));

   //std::cout << split.size() << std::endl;

   std::map<char,int> prob_table;
   for(size_t i = 0;i < split.size();i += 2){
      prob_table.insert(std::make_pair(split[i].at(0),std::stoi(split[i + 1])));
   }

   std::vector<HuffmanNode*> huffman_array;
   CreateHuffmanTree(prob_table,huffman_array);

   std::vector<HuffmanNode> huffman_table;
   CreateHuffmanCode(huffman_table,huffman_array[huffman_array.size()-1],std::string(""));

   bool null_flag = false;
   HuffmanNode *root = huffman_array[huffman_array.size()-1];
   HuffmanNode *it = root;
   while(huf >> str){
      for(size_t i = 0;i < str.size();++i){
	 char c = str[i];
	 if(c == '0'){
	    if(it->L){
	       it = it->L;
	    }
	    else{
	       null_flag = true;
	    }
	 }
	 else if(c == '1'){
	    if(it->R){
	       it = it->R;
	    }
	    else{
	       null_flag = true;
	    }
	 }

	 if(null_flag){
	    rld << it->data;
	    null_flag = false;
	    it = root;
	    --i;
	 }
	 else if(i == str.size() - 1){
	    rld << it->data;
	    null_flag = false;
	    it = root;
	 }
      }
      rld << std::endl;
   }

   huf.close();
}

/*ac.cpp*/
const int ADAPT = 1;
const int NSYM1 = 256;
int loop_num = 0;
void Compression::AriVLC()
{
   ac_encoder ace;
   ac_model acm;

   ac_encoder_init (&ace, (path + std::string(".ari")).c_str());
   ac_model_init (&acm, NSYM1, NULL, ADAPT);

   std::fstream rld(path + std::string(".rlc"),std::ios::in);
   std::string str;
   while(rld >> str){
      for(size_t i = 0;i < str.size();++i){
	 ac_encode_symbol (&ace, &acm, (int)str[i]);
	 ++loop_num;
      }
   }

   ac_encoder_done (&ace);
   ac_model_done (&acm);

   ari_total_byte = (double)ac_encoder_bits (&ace) / 8.0f;
}

void Compression::AriVLD()
{
   ac_decoder acd;
   ac_model acm;
   
   ac_decoder_init (&acd,(path + std::string(".ari")).c_str());
   ac_model_init (&acm, NSYM1, NULL, ADAPT);
   
   std::fstream rlc(path + std::string(".rlc"),std::ios::out);
   while(loop_num){
      char c = (char)ac_decode_symbol (&acd,&acm);
      if(c == 'B'){
	 rlc << c;
	 rlc << std::endl;
      }
      else{
	 rlc << c;
      }
      --loop_num;
   }
   ac_decoder_done (&acd);
   ac_model_done (&acm);
}
