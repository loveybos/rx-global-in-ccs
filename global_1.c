#include <stdio.h>
#include <math.h>
//#include <time.h>
//#include "dsplib.h"
#include "hydicedat.h"

#define SIZE 2768896
#define NUM_ROWS  64
#define NUM_COLS  64
#define BODUAN  169

#pragma SET_DATA_SECTION(".MY_MEM")
float dat_float[SIZE / 4];  				/*float数据*/
float X[BODUAN][NUM_ROWS * NUM_COLS];		/*二维数组*/
//float X_tran[NUM_ROWS * NUM_COLS][BODUAN];  //X'
float rows_sum[BODUAN]; 					/*各波段均值*/
float sigma[BODUAN][BODUAN]; 				/*X的协方差矩阵*/
float d2[NUM_ROWS * NUM_COLS]; 				/*rx算法得出的算子*/
float inv_sigma[BODUAN][BODUAN];			/*sigma的逆矩阵*/
#pragma SET_DATA_SECTION()

void func_read(const unsigned short *p);
void func_to_cov(void);
void func_inv(void);
void func_swap(float*, float*);
void func_to_rx(void);

int main(void)
{
	unsigned int i;
	func_read(dat_hex_h);
	func_to_cov();
	func_inv();
	func_to_rx();
	for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
		printf("%f,", d2[i]);
	printf("end\n");
	return 0;
}

void func_read(const unsigned short *p)
{
	unsigned int i, x;

	printf("正在读取图像...\n");

	// 读取图像数组
	for ( i = 0; i < SIZE; i++ )
	{
		x = ((*p) << 24) | ((*(p + 1)) << 16) | ((*(p + 2)) << 8) | (*(p + 3));
		p += 4;
		dat_float[i] = *(float*)&x;
	}

	printf("图像读取完成!\n");
}


void func_to_cov()
{
	unsigned int i, j, k, flag = 0;
//	clock_t start_t, end_t;
//		   double total_t;
	/*变换成二维数组*/
	printf("正在转化二维数组...\n");
	for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
		for (j = 0; j < BODUAN; j++)
		{
			X[j][i] = dat_float[flag++];
		}
	printf("二维数组转化成功！\n");

	for (i = 0; i < BODUAN; i++)
		rows_sum[i] = 0 ;

	/*求各波段值的和*/
	printf("正在减去均值...\n");
	for (i = 0; i < BODUAN; i++)
		for (j = 0; j < NUM_ROWS * NUM_COLS; j++)
			rows_sum[i] += X[i][j] ;
	/*每个值减去它波段的均值*/
	for (i = 0; i < BODUAN; i++)
		for (j = 0; j < NUM_ROWS * NUM_COLS; j++)
			X[i][j] -= (rows_sum[i] / (NUM_ROWS * NUM_COLS));
	printf("减去均值完成！\n");

	for (i = 0; i < BODUAN; i++)
		for (j = 0; j < BODUAN; j++)
			sigma[i][j] = 0;
	/*求协方差*/
//	 start_t = clock();
//	 DSP_mat_trans(X,BODUAN,NUM_ROWS * NUM_COLS,X_tran);
	printf("正在求协方差...\n");

	for (i = 0; i < BODUAN; i++)
		for (j = 0; j < BODUAN; j++)
		{
			for (k = 0; k < NUM_ROWS * NUM_COLS; k++)
				sigma[i][j] += ((X[i][k]) * (X[j][k]));
			sigma[i][j] = sigma[i][j] / (NUM_ROWS * NUM_COLS);
		}
	printf("求协方差完成！\n");
//	end_t = clock();
//	total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
//	   printf("CPU 占用的总时间：%f\n", total_t  );
}

void func_inv(void)
{

	unsigned int i, j, k;
	float temp, fmax;
	int is[BODUAN] = {0};
	int js[BODUAN] = {0};
	printf("正在求逆...\n");
	for (i = 0; i < BODUAN; i++)
		for (j = 0; j < BODUAN; j++)
			inv_sigma[i][j] = sigma[i][j];

	for (k = 0; k < BODUAN; k++)
	{
		fmax = 0.0;
		for (i = k; i < BODUAN; i++)
		{
			for (j = k; j < BODUAN; j++)
			{
				temp = fabs(inv_sigma[i][j]); /*找最大值*/
				if (temp > fmax)
				{
					fmax = temp;
					is[k] = i; js[k] = j;
				}
			}
		}
		/*if((fmax+1.0)==1.0)
		{
		    return 0;
		}*///////
		if ((i = is[k]) != k)
		{
			for (j = 0; j < BODUAN; j++)
				func_swap(&inv_sigma[k][j], &inv_sigma[i][j]); /*第k行和第i列交换*/
		}
		if ((j = js[k]) != k)
		{
			for (i = 0; i < BODUAN; i++)
				func_swap(&inv_sigma[i][k], &inv_sigma[i][j]); /*第k列和第j列交换*/
		}
		inv_sigma[k][k] = 1.0 / inv_sigma[k][k]; /*求倒*/
		for (j = 0; j < BODUAN; j++)
		{
			if (j != k)
				inv_sigma[k][j] *= inv_sigma[k][k]; /*第k行的每一个元素都除以主对角线上的元素*/
		}
		for (i = 0; i < BODUAN; i++)
		{
			if (i != k)
			{
				for (j = 0; j < BODUAN; j++)
					if (j != k)
						inv_sigma[i][j] = inv_sigma[i][j] - inv_sigma[i][k] * inv_sigma[k][j];
			}

		}
		for (i = 0; i < BODUAN; i++)
		{
			if (i != k)
				inv_sigma[i][k] *= -inv_sigma[k][k];
		}
	}
	for (k = BODUAN - 1; k >= 0; k--)
	{
		if ((j = js[k]) != k)
		{
			for (i = 0; i < BODUAN; i++)
			{
				func_swap(&inv_sigma[j][i], &inv_sigma[k][i]); /*第j行与第k行交换*/
			}
		}
		if ((i = is[k]) != k)
		{
			for (j = 0; j < BODUAN; j++)
			{
				func_swap(&inv_sigma[j][i], &inv_sigma[j][k]); /*第i列与第k列交换*/
			}
		}
	}
	printf("求逆完成!\n");
}

void func_swap(float *a1, float *b1)
{
	float c1;
	c1 = *a1;
	*a1 = *b1;
	*b1 = c1;
}

void func_to_rx(void)
{
	unsigned int i, j, k, flag;
	float d1[BODUAN];

	for (i = 0 ; i < BODUAN ; i++)
		d1[i] = 0;
	for (i = 0; i < NUM_ROWS * NUM_COLS; i++)
		d2[i] = 0;
	printf("正在求rx算子...\n");
	for (flag = 0; flag < NUM_ROWS * NUM_COLS; flag++)
	{
		for (j = 0; j < BODUAN; j++)
		{
			d1[j] = 0;
			for (k = 0; k < BODUAN; k++)
				d1[j] += (X[k][flag] * inv_sigma[k][j]);
		}
		for (j = 0; j < BODUAN; j++)
			d2[flag] += (d1[j] * X[j][flag]);
	}
	printf("rx算子求解成功！\n");
}
