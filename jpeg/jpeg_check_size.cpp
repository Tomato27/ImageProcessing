#include <bits/stdc++.h>
#include <fstream>
using namespace std;

#include "hahuman.h"
#include "zz.h"
#define H 256
#define W 256
#define REP(i, start, end) for(int i=start; i<end; i++)

const char readFileName[]  = "lenna_uchar_256-256.raw";
const char writeFileName[] = "jpeg7.raw";
void readRawFile (const char fname[], const size_t size, const size_t num, void* image);
void copyUchar2Float(unsigned char* src, float* dst, const size_t num);
void writeRawFile(const char fname[], const size_t size, const size_t num, void* image);

//符号器
void DCT(float* f, float* g);
void Quantizer(float* g);
void Encopy_encoder(float* g);
//複合器
void Encopy_decoder(vector<int>&code_ans,float* h);
void Dequantizer(float* h);
void IDCT(float* h,float* o);

float* tmp;

vector<int>code;
vector<int>code_ans;
unsigned char uc[62949/8+1];

int binary(int bina){
  vector<int>rev;
  if(bina>=0){
    for (int i = 0; bina>0 ; i++)
    {
      rev.push_back((bina%2));
      bina = bina/2;
    }
    reverse_copy(rev.begin(), rev.end(), back_inserter(code) );
    rev.clear();
  }
  else{//負のとき
    for (int i = 0; abs(bina)>0 ; i++)
    {
      bina=abs(bina);
      rev.push_back((bina%2)?0:1);
      bina = bina/2;
    }
    reverse_copy(rev.begin(), rev.end(), back_inserter(code) );
    rev.clear();
  }
}

int countx=0,county=0;
int check_q=1;

int main()
{
  // 画像を最初に読み込む領域を確保
  // 読み込む画像がunsigned char型だからunsigned char型の配列を用意する
  // ×float型で計算したいからいきなりfloat型に入れようとするのは間違い
  // 理由:unsigned char型は8bitであるのに対して, float型は32bitだから

  // 画像の読み込み
  unsigned char* image = (unsigned char*)calloc(H * W, sizeof(unsigned char));
  readRawFile(readFileName, sizeof(unsigned char), H * W, image);

  // 画像の配列fにimageをコピー
  // imageはいらなくなるのでfreeしておく
  float* f = (float*)calloc(H * W, sizeof(float));
  copyUchar2Float(image, f, H * W);

  // 回転後の画像を格納する配列gを準備する
  float* g = (float*)calloc(H * W, sizeof(float));

  float* f8 = (float*)calloc(8 * 8, sizeof(float));
  float* g8 = (float*)calloc(8 * 8, sizeof(float));
  tmp = (float*)calloc(8 * 8, sizeof(float));
  FILE	*fpo;	/** ファイルポインタ **/

/************************************************************************************************/
  //符号器
  DCT(f, g);
  Quantizer(g);
  Encopy_encoder(g);

  cout << code.size()<<endl;

  ofstream outputfile("testcode.txt");
  REP(a,0,code.size()){
    int num=code[a];
    outputfile<<num;
  }
  outputfile.close();

  int sum=0;
  int c=0;
  REP(i,0,code.size()){
    REP(j,0,8){
      if(code.size()>i){
        sum+=code[i+j]*pow(2,7-j);
        i++;
      }
    }
    //code_uc.push_back((unsigned char)sum);
    uc[c]=sum;
    //cout<<(int)uc[c];
    c++;
    sum=0;i--;
  }
  //配列にする

  //cout<<code_uc.size();

  ofstream outputfile2("testuc.txt");
  REP(b,0,code.size()/8+1){
    unsigned char c=uc[b];
    outputfile2<<c;
  }
  // REP(b,0,code.size()/4+1){
  //   unsigned char c=code_uc[b];
  //   outputfile2<<c;
  // }
  outputfile2.close();

  float* h = (float*)calloc(H * W, sizeof(float));
  float* o = (float*)calloc(H * W, sizeof(float));
  //複合器

  Encopy_decoder(code_ans,h);
  Dequantizer(h);
  IDCT(h,o);

/************************************************************************************************/
  // writeRawFile(writeFileName, sizeof(float), H * W, o);
  // fclose(fpo);
  //
  // free(image);
  // free(f);
  // free(g);
}

void copyUchar2Float(unsigned char* src, float* dst, const size_t num)
{
    for (size_t i = 0; i < num; i++) { dst[i] = (float)src[i]; }
}

void readRawFile(const char fname[], const size_t size, const size_t num, void* image)
{
    // ファイルを開く
    FILE* fp = fopen(fname, "rb");

    // ファイルを開くことができなかった場合のエラー処理
    if ( NULL == fp )
    {
        printf("failed to open %s\n", fname);
        exit(-1);
    }

    // データの読み込み
    size_t ret = fread(image, size, num, fp);

    // データを読み込むことができなかった場合のエラー処理
    if ( num != ret )
    {
        printf("failed to read %s\n", fname);
        fclose(fp);
        exit(-1);
    }

    // ファイルを閉じる
    fclose(fp);
}

void writeRawFile(const char fname[], const size_t size, const size_t num, void* image)
{
    // ファイルを開く
    FILE* fp = fopen(fname, "wb");

    // ファイルを開くことができなかった場合のエラー処理
    if ( NULL == fp )
    {
        printf("failed to open %s\n", fname);
        exit(-1);
    }

    // データの書き出し
    size_t ret = fwrite(image, size, num, fp);

    // データを書き込むことができなかった場合のエラー処理
    if ( num != ret )
    {
        printf("failed to write %s\n", fname);
        fclose(fp);
        exit(-1);
    }

    // ファイルを閉じる
    fclose(fp);
}

void DCT(float* f, float*g){
  float C1,C2;
  REP(a,0,32){
    REP(b,0,32){
      REP(u,0,8){
        REP(v,0,8){
          float sum=0.0;
          REP(j,0,8){
            REP(k,0,8){
              sum+=f[(j+8*a)*H+k+8*b]*cos((2*j+1)*u*M_PI/16)*cos((2*k+1)*v*M_PI/16);
            }
          }
          if(u==0){C1=1.0/(sqrt(2.0));}
          else{C1=1.0;}
          if(v==0){C2=1.0/(sqrt(2.0));}
          else{C2=1.0;}
          g[(a*8+u)*H+v+8*b]=sum*(C1*C2/4);
        }
      }
    }
  }
}

void Quantizer(float* g){
  REP(a,0,32){
    REP(b,0,32){
      REP(x,0,8){//それぞれの8*8走査
        REP(y,0,8){
          if(g[(x+a*8)*H+y+8*b]>0){
            tmp[x*8+y]=(int)(0.5+g[(x+a*8)*H+y+8*b]*(1/check_q)/q_table[x][y]);
          }
          else{
            tmp[x*8+y]=(int)(g[(x+a*8)*H+y+8*b]*(1/check_q)/q_table[x][y]-0.5);
          }
        }
      }

      //Zigzag scan
      int k=0;
      REP(x1,0,8){
        REP(y1,0,8){
            g[(x1+a*8)*H+y1+8*b]=tmp[zigzag[k]];
            tmp[zigzag[k]]=0;k++;
        }
      }
      if(a==0&&b==25){
        REP(c,0,8){
          REP(d,0,8){
              //cout<<g[(c+a*8)*H+d+8*b]<<" ";
          }
        }
      }
    }
  }
}

void Encopy_encoder(float* g){
  int diff[32*32];//DC成分の誤差
  int run=0;
  REP(a,0,32){
    REP(b,0,32){
      run=0;
      REP(x,0,8){//それぞれの8*8走査
        REP(y,0,8){

          //DC成分
          if(x==0&&y==0){
            if(a==0&&b==0){
              diff[a*32+b]=g[0];
            }
            else{
              //これ左端から右端もOK?
              if(b==0){
                diff[a*32+b]=g[(a*8*H+8*b)]-g[((a-1)*8*H+255-7)];//8*31-1
              }
              else{
                diff[a*32+b]=g[(a*8*H+8*b)]-g[(a*8*H+8*(b-1))];
              }
            }

            //caliculate SSSS
            int ssss=0;
            if(diff[a*32+b]!=0){
              ssss=1+log2(abs(diff[a*32+b]));
              //配列でforのmaxをdc_length_tableにすればよいのでは
              REP(i,0,dc_length_table[ssss]){
                //SSSSの2進表現
                code.push_back(dc_code_table[ssss][i]);
              }

              //diffの2進表現要↓
              if(diff[a*32+b]!=0){
                binary(diff[a*32+b]);
              }
            }
            else{
              //配列でforのmaxをdc_length_tableにすればよいのでは
              REP(i,0,dc_length_table[ssss]){
                //SSSSの2進表現
                code.push_back(dc_code_table[ssss][i]);
              }
            }
          }

          //AC成分
          else{
            int ssss=0;
            if(g[(x+a*8)*H+y+8*b]!=0){
              ssss=1+log2(abs(g[(x+a*8)*H+y+8*b]));
              while(run>15){
                run-=15;
                for(int i1=0;i1<11;i1++){//特殊符号
                  code.push_back(ac_code_table[165][i1]);
                }
                //cout<<code.size()<<" ";
              }
              for(int i2=0;i2<ac_length_table[run*11+ssss];i2++){////11?
                code.push_back(ac_code_table[run*11+ssss][i2]);//////11?
              }
              binary(g[(x+a*8)*H+y+8*b]);
              run=0;
            }
            else if(g[(x+a*8)*H+y+8*b]==0){
              run++;
            }
            if(x==7&&y==7){
              //0のままEOBまで行ったら?→1010のみでOK
              for(int i3=0;i3<4;i3++){
                code.push_back(ac_code_table[0][i3]);
              }
              run=0;
            }
          }
        }
      }
      }
    }
  }

void Encopy_decoder(vector<int> &code,float* h){
  int x=0;
  while(x<code.size()/8+1){
    REP(y,0,8){
      if(uc[x]>=pow(2,7-y)){
        uc[x]-=pow(2,7-y);
        code_ans.push_back(1);
      }
      else{code_ans.push_back(0);}
    }
    x++;
  }
  while(code_ans.size()>code.size()){
    code_ans.pop_back();
  }
  //cout<<"ok";
  int countx=0,county=0;
  int countc=0;
  int a=0;
  int pre_DC=0;
  while(a<=code_ans.size()-1){
    countc++;
    float* in=(float*)calloc(8*8, sizeof(float));
    int out_index=0;//ブロックごとのindex:0=63
    int dc_j=0;
    { //DC成分の処理
        bool check_bits=true;
        REP(i,0,12){/////SSSSの確定
          int b=0;//途中変数
          REP(j,0,dc_length_table[i]){
            check_bits=true;
            if(code_ans[a+j]==dc_code_table[i][j]){//16桁の場合も考慮せねば
              b++;
              if(j==dc_length_table[i]-1){
                a+=dc_length_table[i]-1;////?
                dc_j=j;
                check_bits=true;
                break;
              }
            }
            else{
              if(dc_code_table[i][j]!=-1){
                check_bits=false;
                break;
              }
            }
          }
          if(check_bits==true){
            //追加コードのcheck
            if(code_ans[a+1]==1){//正の値 j,i?
              REP(k,1,i+1){
                in[0]+=pow(2,i-k)*code_ans[a+k];
              }
            }
            else{//負の値
              REP(k,1,i+1){
                in[0]+=-1*pow(2,i-k)*(code_ans[a+k]^1);
              }
            }
            a+=i+1;//+1必要かも
            out_index++;//出力配列のindex更新
            in[0]+=pre_DC;
            break;
          }
        }
    }
    {   //AC成分
        bool check_bits=true;
        int b;
        bool flag_ac=true;
        while(flag_ac){//これ入れる！/////////////////////////////
          REP(i,0,176){
            check_bits=true;
            if(i!=0&&i!=165&&i%11==0){i++;}
            REP(j,0,ac_length_table[i]){//ここもDC同様書き換える
              if(code_ans[a+j]!=ac_code_table[i][j]){
                check_bits=false;
                break; //必ずしも無理してbreakする必要もない?
              }
              else{
                if(j==ac_length_table[i]-1){
                  check_bits=true;
                  break;
              }
            }
          }

            if(check_bits){//ACの何かに合致したら
              a+=ac_length_table[i];
              if(i==0){/////////↓1010の処理では？
                //0で埋める必要あり
                while(out_index<64){
                  in[out_index]=0;
                    out_index++;
                }
                flag_ac=false;
                break;
              }
              else if(i==165){//特殊符号
                REP(z,0,15){
                  in[out_index]=0;
                  out_index++;
                }
                break;
              }
              else{
                REP(x,0,i/11){////////////?run(0の数)
                  in[out_index]=0;
                  out_index++;
                }
                //追加コードのcheck
                if(code_ans[a]==1){
                  REP(k,1,i%11+1){///size(SSSS) a
                    in[out_index]+=pow(2,i%11-k)*code_ans[a];
                    a++;
                  }
                  out_index++;
                  break;
                }
                else{
                  REP(k,1,i%11+1){///run、sizeの後の数字を代入
                    in[out_index]+=-1*pow(2,i%11-k)*(code_ans[a]^1);
                    a++;
                  }
                  out_index++;
                  break;
                }
                b=i;
                break;
              }
            }
          }
        }

        REP(i1,0,8){
          REP(j1,0,8){
            h[(i1+county*8)*H+j1+8*countx]=in[i1*8+j1];
          }
        }
        if(countx==31){
          countx=0;
          county++;
        }
        else if(county>31){
          cout<<county;
        }
        else{
          countx++;
        }
        pre_DC=in[0];
        free(in);
    }
 }
}

void Dequantizer(float* h){
  REP(a,0,32){
    REP(b,0,32){
      //inv_Zigzag scan
      int k=0;
      REP(x1,0,8){
        REP(y1,0,8){
          tmp[((int)(zigzag[k]/8)*8)+(int)(zigzag[k]%8)]=h[(a*8+x1)*H+b*8+y1];
          k++;
        }
      }
      REP(x,0,8){//それぞれの8*8走査
        REP(y,0,8){
          if(tmp[x*8+y]>0){
            h[(a*8+x)*H+b*8+y]=(int)(check_q*tmp[x*8+y]*q_table[x][y]-0.5);
          }
          else{
            h[(a*8+x)*H+b*8+y]=(int)(check_q*tmp[x*8+y]*q_table[x][y]+0.5);
          }
        }
      }
    }
  }
}

void IDCT(float* h,float* o){
  float C1,C2;
  REP(a,0,32){
    REP(b,0,32){
      REP(j,0,8){
        REP(k,0,8){
            float sum=0.0;
            REP(u,0,8){
                REP(v,0,8){
                  if(u==0){C1=1.0/(sqrt(2.0));}
                  else{C1=1.0;}
                  if(v==0){C2=1.0/(sqrt(2.0));}
                  else{C2=1.0;}
                  sum+=C1*C2*h[(u+8*a)*H+v+8*b]*cos((2*j+1)*u*M_PI/16)*cos((2*k+1)*v*M_PI/16);
                }
            }
            o[(a*8+j)*H+k+8*b]=sum/4;
        }
      }
    }
  }
}
