/*
 *                        _oo0oo_
 *                       o8888888o
 *                       88" . "88
 *                       (| -_- |)
 *                       0\  =  /0
 *                     ___/`---'\___
 *                   .' \\|     |// '.
 *                  / \\|||  :  |||// \
 *                 / _||||| -:- |||||- \
 *                |   | \\\  - /// |   |
 *                | \_|  ''\---/''  |_/ |
 *                \  .-\__  '-'  ___/-. /
 *              ___'. .'  /--.--\  `. .'___
 *           ."" '<  `.___\_<|>_/___.' >' "".
 *          | | :  `- \`.;`\ _ /`;.`/ - ` : | |
 *          \  \ `_.   \_ __\ /__ _/   .-` /  /
 *      =====`-.____`.___ \_____/___.-`___.-'=====
 *                        `=---='
 * 
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */

/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/



#include "boundaryconditions.hpp"
#include "ghost.hpp"
#include <cmath>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__device__ void calculateSymmetricPointLine(float A, float B, int ki, int j, int n, float &xSym, float &ySym) {
    // 给定点的坐标
    float x = ki;
    float y = j + n + 1;
    
    // 计算垂足坐标 x0 和 y0
    float x0 = (x + A * y - A * B) / (A * A + 1);
    float y0 = (A * A * y + A * x + B) / (A * A + 1);
    
    // 计算对称点的坐标 xSym 和 ySym
    xSym = 2 * x0 - x;
    ySym = 2 * y0 - y;
}

__device__ void calculateSymmetricPointCircle(float x_o, float y_o, int ki, int kj, float R, float &xSym, float &ySym) {
    // 给定点的坐标
    float x_p = ki;
    float y_p = kj;
    float dx = x_p - x_o;
    float dy = y_p - y_o;
	float d = sqrt(dx*dx + dy*dy);
	xSym = dx * (2 * R - d)/d + x_o;
    ySym = dy * (2 * R - d)/d + y_o;
}

__device__ void calculateDistance(float xSym, float ySym, int x, int y, float &distance) {
	distance = sqrt(pow((xSym - x),2) + pow((ySym - y),2));
}
	
__device__ void calculatePerpendicularAngle(float xSym, float ySym, float x_o, float y_o, float &rad_o) {
    float slope = (ySym - y_o) / (xSym - x_o);
	float tan_slope = -1/slope;
    // 计算OP与X轴的夹角θ
    float theta = atan(tan_slope);
    // 计算垂线与X轴的夹角
    rad_o = -theta;
}	

template<BoundaryConditions BCs, bool XDIR, bool downstream>
__global__ void setBoundaryConditionsKernel(Mesh<GPU>::type u) {
  //const bool YDIR = !XDIR, upstream = !downstream;

  const int i = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();

  if (XDIR && u.exists(i, 0)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
         //u(i, u.activeNy() + 1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, u.activeNy() - 2, k);
        //u(i, u.activeNy()    , k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, u.activeNy() - 1, k);
		//u(i, -2, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 1, k);
        //u(i, -1, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 0, k); 
      // u(i, u.activeNy() + 1, k) = u(i, 1, k);
      // u(i, u.activeNy()    , k) =  u(i, 0, k);
       u(i, u.activeNy() + 1, k) = u(i, u.activeNy() - 2, k); 
       u(i, u.activeNy()    , k) = u(i, u.activeNy() - 1, k);

      //} else {
	  }
  }
  }
	  if (XDIR && u.exists(i, 0)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
		  //u(i, -2, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 1, k);
         //u(i, -1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 0, k);
		//u(i, u.activeNy() + 1, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, u.activeNy() - 2, k);
        //u(i, u.activeNy()    , k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, u.activeNy() - 1, k);
        //u(i, -2, k) = u(i, u.activeNy() - 2, k); 
        //u(i, -1, k) = u(i, u.activeNy() - 1, k);
       u(i, -2, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 1, k);
       u(i, -1, k) = u(i, 0, k);

      }
    }
	  // u(i + wedge_x * cell_x , u.activeNy() + 1, YMOMENTUM) = - u(i + wedge_x * cell_x, u.activeNy() -2, YMOMENTUM); 
      // u(i + wedge_x * cell_x, u.activeNy() , YMOMENTUM) = - u(i + wedge_x * cell_x, u.activeNy() -1, YMOMENTUM);
  }
  if (!XDIR && u.exists(0, i)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
         //u(-2, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(1, i, k);
         //u(-1, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(0, i, k);
        u(-2, i, k) = u(1, i, k);
        u(-1, i, k) = u(0, i, k); 
      } 
  }}
	  if (!XDIR && u.exists(0, i)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
        // u(u.activeNx() + 1, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(u.activeNx() - 2, i, k);
        //u(u.activeNx()    , i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(u.activeNx() - 1, i, k);
        u(u.activeNx() + 1, i, k) = u(u.activeNx() - 2, i, k);
       u(u.activeNx()    , i, k) = u(u.activeNx() - 1, i, k);
      }
    }
  }


const int ki = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells(); 
const int kj = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();

float rad_theta_C = theta_C * PI / 180; 
float tan_theta_C = tan(rad_theta_C);
float cos_theta_C = cos(rad_theta_C);
float sin_theta_C = sin(rad_theta_C);

//Semi-cylindrical bump surface////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	float x_o = (l_bump * cos_theta_C + A_x) * cell_x;
	float y_o = (AB_y - sin_theta_C * l_bump) * cell_y;

	float R = r_bump * cell_x;
	float xSym,ySym,D11,D12,D21,D22,x1,x2,y1,y2, rad_o;
	
	if (XDIR && u.exists(ki, 0) && ki >= (x_o - R/1.414) && ki <= (x_o + R/1.414)) {
		for (int n = 0; n < 2; n++) {
		for (int j = y_o - R - 1; j < y_o + R + 1; j++) {
		if (((pow((ki - x_o),2) + pow((j - y_o),2)) >= pow(R,2)) && ((pow((ki - x_o),2) + pow((j - y_o),2)) <= pow(R + 1.415,2))){
			//if (j <= (y_o - R/1.414)){	//top
			if (j >= (y_o)){
			    int kj = j + n;
				if ((pow((ki - x_o),2) + pow((kj - y_o),2)) >= pow(R,2) && kj >= y_o){
				calculateSymmetricPointCircle(x_o, y_o, ki, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				          if (D11 < 0.00001) {
            D11 = 0.00001;
          }
          if (D12 < 0.00001) {
            D12 = 0.00001;
          }
          if (D21 < 0.00001) {
            D21 = 0.00001;
          }
          if (D22 < 0.00001) {
            D22 = 0.00001;
          }
				float Dsum = D11 + D12 + D21 + D22; 
				
			     for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				calculatePerpendicularAngle(xSym, ySym, x_o, y_o, rad_o);
				
				float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
/*
				if ((pow((x1 - x_o),2) + pow((y2 - y_o),2)) <= pow(R,2)){
                u(ki, kj, DENSITY) = 100;
				}
*/
					
				
	           	float v_tan = vx_temp * cos(rad_o) - vy_temp * sin(rad_o);
	            float v_norm = vx_temp * sin(rad_o) + vy_temp * cos(rad_o);	
	            //float v_x_temp = -v_norm * sin(rad_o) + v_tan * cos(rad_o); //slip wall
	            //float v_y_temp = -v_norm * cos(rad_o) - v_tan * sin(rad_o); //slip wall
				float v_x_temp = -v_norm * sin(rad_o) - v_tan * cos(rad_o); 
	            float v_y_temp = -v_norm * cos(rad_o) + v_tan * sin(rad_o); 
	  
	            u(ki, kj, XMOMENTUM) = v_x_temp * u(ki, kj, DENSITY);
	            u(ki, kj, YMOMENTUM) = v_y_temp * u(ki, kj, DENSITY); 
				
				 
			}		
			}
			//else if (j >= (y_o + R/1.414)){	//bot
			else if (j <= (y_o)){	//bot
			//if (j <= (y_o)){
			    int kj = j - n;
				if ((pow((ki - x_o),2) + pow((kj - y_o),2)) >= pow(R,2) && kj < y_o){
				calculateSymmetricPointCircle(x_o, y_o, ki, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				          if (D11 < 0.00001) {
            D11 = 0.00001;
          }
          if (D12 < 0.00001) {
            D12 = 0.00001;
          }
          if (D21 < 0.00001) {
            D21 = 0.00001;
          }
          if (D22 < 0.00001) {
            D22 = 0.00001;
          }
				float Dsum = D11 + D12 + D21 + D22; 
				//u(ki, kj, DENSITY) = 10;
	           // u(ki, kj, DENSITY) = 10; 
				
			     for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				calculatePerpendicularAngle(xSym, ySym, x_o, y_o, rad_o);
	            float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
	           	float v_tan = vx_temp * cos(rad_o) - vy_temp * sin(rad_o);
	            float v_norm = vx_temp * sin(rad_o) + vy_temp * cos(rad_o);	
	            float v_x_temp = -v_norm * sin(rad_o) - v_tan * cos(rad_o); 
	            float v_y_temp = -v_norm * cos(rad_o) + v_tan * sin(rad_o); 
				//float v_x_temp = -v_norm * sin(rad_o) + v_tan * cos(rad_o);// slip wall
	            //float v_y_temp = -v_norm * cos(rad_o) - v_tan * sin(rad_o);// slip wall
	  
	            u(ki, kj, XMOMENTUM) = v_x_temp * u(ki, kj, DENSITY);
	            u(ki, kj, YMOMENTUM) = v_y_temp * u(ki, kj, DENSITY); 
				 
			}		
			}
			/*
			else if (ki <= (x_o)){
			    int kj = j;
				int ki_temp = ki + n + 1;
				calculateSymmetricPointCircle(x_o, y_o, ki_temp, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				float Dsum = D11 + D12 + D21 + D22; 
				//u(ki, kj, DENSITY) = 10;
	           // u(ki, kj, DENSITY) = 10; 
				
			     for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki_temp, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				calculatePerpendicularAngle(xSym, ySym, x_o, y_o, rad_o);
	            float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
	           	float v_tan = vx_temp * cos(rad_o) - vy_temp * sin(rad_o);
	            float v_norm = vx_temp * sin(rad_o) + vy_temp * cos(rad_o);	
	            float v_x_temp = -v_norm * sin(rad_o) + v_tan * cos(rad_o); 
	            float v_y_temp = -v_norm * cos(rad_o) - v_tan * sin(rad_o); 
	  
	            u(ki_temp, kj, XMOMENTUM) = v_x_temp * u(ki_temp, kj, DENSITY);
	            u(ki_temp, kj, YMOMENTUM) = v_y_temp * u(ki_temp, kj, DENSITY);
			}
						else if (ki >= (x_o)){
			    int kj = j;
				int ki_temp = ki - n;
				calculateSymmetricPointCircle(x_o, y_o, ki_temp, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				float Dsum = D11 + D12 + D21 + D22; 
				//u(ki, kj, DENSITY) = 10;
	           // u(ki, kj, DENSITY) = 10; 
				
			     for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki_temp, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				calculatePerpendicularAngle(xSym, ySym, x_o, y_o, rad_o);
	            float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
	           	float v_tan = vx_temp * cos(rad_o) - vy_temp * sin(rad_o);
	            float v_norm = vx_temp * sin(rad_o) + vy_temp * cos(rad_o);	
	            float v_x_temp = -v_norm * sin(rad_o) + v_tan * cos(rad_o); 
	            float v_y_temp = -v_norm * cos(rad_o) - v_tan * sin(rad_o); 
	  
	            u(ki_temp, kj, XMOMENTUM) = v_x_temp * u(ki_temp, kj, DENSITY);
	            u(ki_temp, kj, YMOMENTUM) = v_y_temp * u(ki_temp, kj, DENSITY);
			}
			*/
			/*
            else if (kj >= (y_o + R/1.414)){
                int j = kj - n;			
				calculateSymmetricPointCircle(x_o, y_o, ki, j, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				float Dsum = D11 + D12 + D21 + D22; 
			    for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki, j, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				
	            float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
	           	float v_tan = vx_temp * cos(rad_theta_C) - vy_temp * sin(rad_theta_C);
	            float v_norm = vx_temp * sin(rad_theta_C) + vy_temp * cos(rad_theta_C);	
	            float v_x_temp = -v_norm * sin(rad_theta_C) + v_tan * cos(rad_theta_C); 
	            float v_y_temp = -v_norm * cos(rad_theta_C) - v_tan * sin(rad_theta_C); 
	  
	            u(ki, j, XMOMENTUM) = v_x_temp * u(ki, j, DENSITY);
	            u(ki, j, YMOMENTUM) = v_y_temp * u(ki, j, DENSITY); 
			}
			*/
			/*
			else if (ki <= (x_o) && kj <= (y_o - R/1.414)){
				int i = ki + n + 1;
				calculateSymmetricPointCircle(x_o, y_o, i, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				float Dsum = D11 + D12 + D21 + D22; 
			    for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(i, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				
	            float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
	           	float v_tan = vx_temp * cos(rad_theta_C) - vy_temp * sin(rad_theta_C);
	            float v_norm = vx_temp * sin(rad_theta_C) + vy_temp * cos(rad_theta_C);	
	            float v_x_temp = -v_norm * sin(rad_theta_C) + v_tan * cos(rad_theta_C); 
	            float v_y_temp = -v_norm * cos(rad_theta_C) - v_tan * sin(rad_theta_C); 
	  
	            u(i, kj, XMOMENTUM) = v_x_temp * u(i, kj, DENSITY);
	            u(i, kj, YMOMENTUM) = v_y_temp * u(i, kj, DENSITY); 

			}
			*/
			/*
			if (kj <= (y_o - R/1.414)){		
				u(x1, y1, DENSITY) = 2;
				u(x2, y2, DENSITY) = 2;
				u(x1, y2, DENSITY) = 2;
				u(x2, y1, DENSITY) = 2;
				//u(ki, kj, DENSITY) = 2;
			}
			else if (kj >= (y_o + R/1.414)){
				u(x1, y1, DENSITY) = 4;
				u(x2, y2, DENSITY) = 4;
				u(x1, y2, DENSITY) = 4;
				u(x2, y1, DENSITY) = 4;
			}
			else if (ki <= (x_o)){
				u(x1, y1, DENSITY) = 6;
				u(x2, y2, DENSITY) = 6;
				u(x1, y2, DENSITY) = 6;
				u(x2, y1, DENSITY) = 6;
			}
			else if (ki >= (x_o)){
				u(x1, y1, DENSITY) = 8;
				u(x2, y2, DENSITY) = 8;
				u(x1, y2, DENSITY) = 8;
				u(x2, y1, DENSITY) = 8;
			}
			else{
				u(ki, kj, DENSITY) = 10;
			}*/
			}
		}
		}
		/*
		if ((R * R - pow((ki - x_o),2))>= 0){
		    float sqrtPart = sqrt(R * R - pow((ki - x_o),2));	
		    float y_upper_half = sqrtPart + y_o;
		    float y_lower_half = - sqrtPart + y_o;
		    u(ki, y_upper_half, DENSITY) = 10;
			u(ki, y_lower_half, DENSITY) = 5;
		}
		*/
	}

	if (!XDIR && u.exists(0, kj) && kj >= (y_o - R/1.414) && kj <= (y_o + R/1.414)) {
		for (int n = 0; n < 2; n++) {
		for (int i = x_o - R - 1; i < x_o + R + 1; i++) {
		if (((pow((i - x_o),2) + pow((kj - y_o),2)) >= pow(R,2)) && ((pow((i - x_o),2) + pow((kj - y_o),2)) <= pow(R + 1.415,2))){
			//if (j <= (y_o - R/1.414)){	//left
			if (i <= (x_o)){
			    int ki = i + n;
				if ((pow((ki - x_o),2) + pow((kj - y_o),2)) >= pow(R,2) && ki <= x_o){
				calculateSymmetricPointCircle(x_o, y_o, ki, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				          if (D11 < 0.00001) {
            D11 = 0.00001;
          }
          if (D12 < 0.00001) {
            D12 = 0.00001;
          }
          if (D21 < 0.00001) {
            D21 = 0.00001;
          }
          if (D22 < 0.00001) {
            D22 = 0.00001;
          }
				float Dsum = D11 + D12 + D21 + D22; 
				
			     for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				calculatePerpendicularAngle(xSym, ySym, x_o, y_o, rad_o);
				
				float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
/*
				if ((pow((x1 - x_o),2) + pow((y2 - y_o),2)) <= pow(R,2)){
                u(ki, kj, DENSITY) = 100;
				}
*/
					
				
	           	float v_tan = vx_temp * cos(rad_o) - vy_temp * sin(rad_o);
	            float v_norm = vx_temp * sin(rad_o) + vy_temp * cos(rad_o);	
	            float v_x_temp = -v_norm * sin(rad_o) + v_tan * cos(rad_o); //slip wall
	            float v_y_temp = -v_norm * cos(rad_o) - v_tan * sin(rad_o); //slip wall
				//float v_x_temp = -v_norm * sin(rad_o) - v_tan * cos(rad_o); 
	            //float v_y_temp = -v_norm * cos(rad_o) + v_tan * sin(rad_o); 
	  
	            u(ki, kj, XMOMENTUM) = v_x_temp * u(ki, kj, DENSITY);
	            u(ki, kj, YMOMENTUM) = v_y_temp * u(ki, kj, DENSITY); 
				
				 
			}		
			}
			else if (i >= (x_o)){	//right
			//if (j <= (y_o)){
			    int ki = i - n;
				if ((pow((ki - x_o),2) + pow((kj - y_o),2)) >= pow(R,2) && ki > x_o){
				calculateSymmetricPointCircle(x_o, y_o, ki, kj, R, xSym, ySym);
	            x1 = floor(xSym); x2 = ceil(xSym); y1 = floor(ySym); y2 = ceil(ySym);	
			    calculateDistance(xSym,ySym,x1,y1,D11);
	            calculateDistance(xSym,ySym,x1,y2,D12);
	            calculateDistance(xSym,ySym,x2,y1,D21);
	            calculateDistance(xSym,ySym,x2,y2,D22);
				          if (D11 < 0.00001) {
            D11 = 0.00001;
          }
          if (D12 < 0.00001) {
            D12 = 0.00001;
          }
          if (D21 < 0.00001) {
            D21 = 0.00001;
          }
          if (D22 < 0.00001) {
            D22 = 0.00001;
          }
				float Dsum = D11 + D12 + D21 + D22; 
				//u(ki, kj, DENSITY) = 10;
	           // u(ki, kj, DENSITY) = 10; 
				
			     for (int k = 0; k < NUMBER_VARIABLES; k++) {
		        u(ki, kj, k) = (D11/Dsum) * u(x1,y1,k) + (D12/Dsum) * u(x1,y2,k) + (D21/Dsum) * u(x2,y1,k) + (D22/Dsum) * u(x2,y2,k);
		
	            }
				calculatePerpendicularAngle(xSym, ySym, x_o, y_o, rad_o);
	            float vx_temp = (D11/Dsum) * u(x1,y1,XMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,XMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,XMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,XMOMENTUM)/u(x2,y2,DENSITY);
	            float vy_temp = (D11/Dsum) * u(x1,y1,YMOMENTUM)/u(x1,y1,DENSITY) + (D12/Dsum) * u(x1,y2,YMOMENTUM)/u(x1,y2,DENSITY) + (D21/Dsum) * u(x2,y1,YMOMENTUM)/u(x2,y1,DENSITY) + (D22/Dsum) * u(x2,y2,YMOMENTUM)/u(x2,y2,DENSITY);
	           	float v_tan = vx_temp * cos(rad_o) - vy_temp * sin(rad_o);
	            float v_norm = vx_temp * sin(rad_o) + vy_temp * cos(rad_o);	
	            //float v_x_temp = -v_norm * sin(rad_o) - v_tan * cos(rad_o); 
	            //float v_y_temp = -v_norm * cos(rad_o) + v_tan * sin(rad_o); 
				float v_x_temp = -v_norm * sin(rad_o) + v_tan * cos(rad_o);// slip wall
	            float v_y_temp = -v_norm * cos(rad_o) - v_tan * sin(rad_o);// slip wall
	  
	            u(ki, kj, XMOMENTUM) = v_x_temp * u(ki, kj, DENSITY);
	            u(ki, kj, YMOMENTUM) = v_y_temp * u(ki, kj, DENSITY); 
				 
			}		
			}
			}
		}
		}
	}
	
  }

