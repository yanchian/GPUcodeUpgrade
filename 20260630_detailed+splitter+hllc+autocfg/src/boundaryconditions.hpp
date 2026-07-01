/*
 * @Author: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 * @Date: 2024-10-25 16:50:08
 * @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 * @LastEditTime: 2024-10-25 16:53:33
 * @FilePath: \undefinedc:\Users\2211\Documents\NetSarang Computer\8\Xftp\Temporary\boundaryconditions.hpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#pragma once

enum Directions {YDIRECTION, XDIRECTION};
enum Streams {UPSTREAM, DOWNSTREAM};

template<BoundaryConditions BCs, bool XDIR, bool downstream>
__global__ void setBoundaryConditionsKernel(Mesh<GPU>::type u);

template<bool XDIR>
__global__ void setSpecialBoundaryConditionsKernel(Mesh<GPU>::type u);

