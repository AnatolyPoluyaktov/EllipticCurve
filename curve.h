#ifndef NULL
#define NULL (void*)0
#endif

#ifndef CURVE_H
#define CURVE_H
#include <openssl/bn.h>

// Набор параметров id-tc26-gost-3410-2012-256-paramSetA:
#define p_str               "115792089237316195423570985008687907853269984665640564039457584007913129639319"

#define a_str               "87789765485885808793369751294406841171614589925193456909855962166505018127157"
#define x_base_str          "65987350182584560790308640619586834712105545126269759365406768962453298326056"
#define y_base_str          "22855189202984962870421402504110399293152235382908105741749987405721320435292"
#define q_str               "28948022309329048855892746252171976963338560298092253442512153408785530358887"

// Значение, вычисленное на основе вышеперечисленных параметров с помощью Wolfram Mathematica
#define theta_str           "454069018412434321972378083527459607666454479745512801572100703902391945898"




// структура для хранения параметров в виде больших чисел
struct param{
    BIGNUM* p;
    BIGNUM* a;
    BIGNUM* x_base;
    BIGNUM* y_base;
    BIGNUM* q;
    BIGNUM* theta;
};
// структура для хранения параметров квадрики и порождающих элементов.
struct JacobiCurve{
    BIGNUM* Y;
    BIGNUM* X;
    BIGNUM* e;
    BIGNUM* d;
    BIGNUM* Z;
    BIGNUM* p;
};
// структура для хранения точек
struct Point{
    BIGNUM* X;
    BIGNUM* Y;
    BIGNUM* Z;
};

//инициализировать структуру с параметрами
void InitParam(struct param* param);

// инициализировать параметры кривой
void InitJacobiCurve(struct JacobiCurve* curve, struct param* param);


// инициалзиция точки
void InitPoint(struct Point* point, char* X, char* Y, char* Z);

//реализация сложения двух точек
void AddPoints(struct Point* P1, struct Point* P2, struct Point* P3, struct JacobiCurve* curve);

// вычисление кратной точки с помощью алгоритма "лесенка Монтгомери"
void CalculateDegree(struct Point* kP, struct Point* P, struct JacobiCurve* curve, BIGNUM* degree);

// перевод в аффинные координаты
void CastPointToAffine(struct Point* affpoint, struct Point* P, struct JacobiCurve* curve);

// вывод в аффинных координатах
void PrintInAffin(struct Point* point, struct JacobiCurve* curve);

// вывод в проективных
void PrintProjective(struct Point* point);

// проверка , что точка лежит на кривой
int CheckPointIsOnCurve (struct Point* P, struct JacobiCurve* curve);

// проверка, что точки равны
int IsEqual(struct Point* P1, struct Point* P2, struct JacobiCurve* curve);

//получение обратной точки
void GetNegativePoint(struct Point *res, struct Point *point);

// очистка памяти
void FreeParam(struct param* para);
void FreeJacobiCurve(struct JacobiCurve* curve);
void FreePoint(struct Point* P);
#endif // CURVE_H
