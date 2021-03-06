// HW4_黃子青_E24076938_陳俊瑋E24076603.cpp: 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include "CImg.h"
#include<iostream>
#include<vector>
#include<fstream>
using namespace cimg_library;
using namespace std;
CImg<unsigned char> src_image("input_binary_image.bmp");//read in an image
int width = src_image.width(); //get width of the image 600
int height = src_image.height(); //get height of the image 600
unsigned char* ptr;//指向圖片的指標
bool start = true;//開始讀黑的判斷
int black_num = 0;//暫存
int record_column = 0;
vector<vector<int>> record(height);//壓縮檔用的
vector<vector<int>> record2(height);//讀壓縮檔用的
vector<vector<int>> grayValue(height);//存灰階的值
int *segment_num = new int [height] {0};//壓縮檔用的
int *segment_num_2 = new int [height] {0};//讀壓縮檔用的
int main()
{



	
	for (int i = 0; i < height; i++)
		grayValue[i].resize(width);
	//壓縮
	for (int i = 0; i < height; i++)
	{
		record_column = 0;
		for (int j = 0; j < width; j++)
		{
			ptr = src_image.data(j, i);
			grayValue[i][j] = *ptr; //255是白色 0是黑色
			if (grayValue[i][j] == 0 && start == true && j != (width - 1)) //開始讀黑
			{
				record[i].push_back(j);
				record_column = j;
				start = false;
				black_num = 0;
			}


			else if (j == (width -1) && grayValue[i][j] == 0&& start==true) { //只讀到一個黑且是最後一個
				record[i].push_back(j + black_num);
				record[i].push_back(j + black_num);
				black_num = 0;
				segment_num[i] += 2;
				start = true;

			}

			else if (j == (width - 1) && grayValue[i][j] == 0) {    //連續讀到最後
				black_num++;
				record[i].push_back(record_column + black_num);
				
				black_num = 0;
				segment_num[i] += 2;
				start = true;

			}

			else if (grayValue[i][j] == 0 && start == false)//連續讀黑
			{
				black_num++;
			}
			else if (grayValue[i][j] == 255 && start == false)//遇到白的，停止讀黑，把column放到record
			{
				record[i].push_back(record_column + black_num);
				black_num = 0;
				segment_num[i] += 2;
				start = true;

			}

		}
		//cout << endl;
	}
	//輸出壓縮檔，分別輸出一個Number of segment 跟startpoint and endpoint
	ofstream fout;//輸出number of segment
	ofstream fout2;//輸出startpoint and endpoint
	fout.open("Number_of_segments_in_each_row.txt");
	fout2.open("compressed_version.txt");
	CImg<unsigned char> out_img(width, height, 1, 1);
	out_img.fill(255);//底是白色
	
	for (int i = 0; i < height; i++)
	{
		fout << segment_num[i] << endl;
		for (int j = 0; j < record[i].size(); j++)
		{
			fout2 << record[i][j] << " ";

		}
		fout2 << endl;
	}
	
	//讀入壓縮檔，並解壓縮
	ifstream inputFile;
	ifstream inputFile2;
	inputFile.open("Number_of_segments_in_each_row.txt");
	inputFile2.open("compressed_version.txt");

	
	
	for (int i = 0; i < height; i++)
	{
		inputFile >> segment_num_2[i];
		
	
			for (int z = 0; z < segment_num_2[i]; z++)
			{
				int temp;
				inputFile2 >> temp;
				record2[i].push_back(temp);
			}
		

			for (int j = 0; j < segment_num_2[i]; j += 2)
			{
				if ((j + 1 < record2[i].size())) //bound checking
					for (int column = record2[i][j]; (column <= record2[i][j + 1]); column++) {
						unsigned char *ptr = out_img.data(column,i);
						*ptr = int(0);//讀到黑色的起始位置跟終點位置=>把經過的column設成黑色
					}

			}
		
	
	}
	
	//save an image with file name:"output_image.bmp"
	out_img.save("output_image.bmp");
	//open a window to display an image on the screen with title : "output"
	CImgDisplay main_disp(out_img, "Output");
	// wait until main window is closed
	while (!main_disp.is_closed()) {
		main_disp.wait();
	}
}