#include <openssl/bn.h>
#include "curve.h"
#include <stdio.h>
#include <string.h>

void Testing(){
    struct param prm = {NULL, NULL, NULL, NULL, NULL, NULL};
    InitParam(&prm);
    //вывожу наборы параметров
    printf(" Набор параметров  id-tc26-gost-3410-2012-256-paramSetA:\n");
    printf("p=%s\n", BN_bn2dec(prm.p));
    printf("a=%s\n", BN_bn2dec(prm.a));
    printf("x_base=%s\n", BN_bn2dec(prm.x_base));
    printf("y_base=%s\n", BN_bn2dec(prm.y_base));
    printf("q=%s\n\n", BN_bn2dec(prm.q));
    printf("посчитали тета параметр:\n");
    printf("theta=%s\n\n", BN_bn2dec(prm.theta));


    //параметры квадрики
    printf("параметры Квадрики Якоби:\n");
    struct JacobiCurve curve={NULL, NULL, NULL, NULL, NULL, NULL};
    InitJacobiCurve(&curve, &prm);
    printf("e=%s\n", BN_bn2dec(curve.e));
    printf("d=%s\n", BN_bn2dec(curve.d));
    printf("X_base=%s\n", BN_bn2dec(curve.X));
    printf("Y_base=%s\n", BN_bn2dec(curve.Y));
    printf("Z_base=%s\n\n", BN_bn2dec(curve.Z));

    //тестирую работу с нейтральным элементом
    printf("ТЕСТ 1:\n");
    struct Point E = {NULL, NULL, NULL};
    InitPoint(&E ,"0", "1", "1");
    printf("Нейтральный элемент Е:\n");
    printf("в проективных координатах:\n");
    PrintProjective(&E);
    printf("в афинных:\n");
    PrintInAffin(&E, &curve);

    if(CheckPointIsOnCurve(&E, &curve)){
        printf("точка E находится на кривой\n\n");
    }
    else{
        printf("точка E не находится на кривой\n\n");
    }


    // тестирую функций с порождающим элементом
    printf("ТЕСТ 2:\n");
    struct Point P_base = {NULL, NULL, NULL};

    InitPoint(&P_base,"0","1","1");
    BN_copy(P_base.X,curve.X);
    BN_copy(P_base.Y,curve.Y);
    BN_copy(P_base.Z,curve.Z);
    printf("порождающий элемент P_base в аффинных координатах:\n");
    PrintInAffin(&P_base,&curve);
    if(CheckPointIsOnCurve(&P_base,&curve)){
        printf("точка P_base находится на кривой\n\n");
    }
    else{
        printf("точка P_base не находится на кривой\n\n");
    }

    // тест E+P=P
    printf("ТECT 3:\n");
    printf("проверим, что E+P_base = P_base\n");
    struct Point P1 = {NULL, NULL, NULL};
    InitPoint(&P1, "2", "2", "2");
    AddPoints(&E, &P_base, &P1, &curve);
    if (!(IsEqual(&P_base, &P1, &curve))){
        printf("E+P и P равны\n\n");
    }
    else{
        printf("E+P и P не равны\n\n");
    }
    //тест выбрать любую точку и проверить лежит ли она кривой
    printf("ТECT 4:\n");
    printf("проверим пренадлежит ли точка P2=(3:6:8) кривой\n");
    printf("P3 в аффинных:\n");
    struct Point P2 = {NULL, NULL, NULL};
    InitPoint(&P2, "3", "6", "8");
    PrintInAffin(&P2, &curve);
    if(CheckPointIsOnCurve(&P2,&curve)){
        printf("точка P2 находится на кривой\n\n");
    }
    else{
        printf("точка P2 не находится на кривой\n\n");
    }

       //тест qP = E
      printf("ТECT 5:\n");
      printf("проверим, что qP = E\n");
      struct Point resPoint = {NULL, NULL, NULL};
      InitPoint(&resPoint, "1", "1" ,"1");
      CalculateDegree(&resPoint, &P_base, &curve, prm.q);
      printf("нейтральный элемент в аффинных координатах:\n");
      PrintInAffin(&E, &curve);
      printf("qP в аффинных координатах:\n");
      PrintInAffin(&resPoint, &curve);

      //тест [q+1]P = P и [q-1] = -P
      printf("ТECT 6:\n");
      printf("проверим, что [q+1]P = P и [q-1] = -P\n");
      BIGNUM* degree = BN_new();
      BN_dec2bn(&degree, "1");// degree=1;
      BN_add (degree, degree, prm.q);//degree=q+1
      CalculateDegree(&resPoint, &P_base, &curve, degree);
      printf("[q+1]P:\n");
      PrintInAffin(&resPoint, &curve);
      printf("P:\n");
      PrintInAffin(&P_base, &curve);
      if(!(IsEqual(&resPoint, &P_base,&curve))){
                printf("[q+1]P равно P\n\n");
       }
      else{
          printf("[q+1]P не равно P\n\n");
      }

      BN_dec2bn (&degree, "1");                      // degree = 1
      BN_sub (degree, prm.q, degree);                // degree = q-1
      CalculateDegree(&resPoint, &P_base, &curve,degree);
      printf("[q-1]P:\n");
      PrintInAffin(&resPoint, &curve);
      printf("-P:\n");
      struct Point negP = {NULL, NULL, NULL};
      InitPoint(&negP, "3", "6", "8");
      GetNegativePoint(&negP, &P_base);
      PrintInAffin(&negP, &curve);
      if(!(IsEqual(&resPoint, &negP, &curve))){
                printf("[q-1]P равно -P\n\n");
       }
      else{
          printf("[q+1]P не равно P\n\n");
      }



     // придумать любое k , такое что 0 <= k <q  и посчитать [k]P и проверить пренадлежность кривой
     printf("ТECT 7:\n");
     printf("вычислим кратную точку [k]P. пусть k = 100,P = P_base\n");
     BN_dec2bn(&degree, "100");
     CalculateDegree(&resPoint, &P_base, &curve,degree);
     PrintInAffin(&resPoint, &curve);
     if(CheckPointIsOnCurve(&resPoint, &curve)){
         printf("точка [k]P находится на кривой\n\n");
     }
     else{
         printf("точка [k]P не находится на кривой\n\n");
     }



       //сгенерировать случайное k , такое что 0 <=k < q и проверить пренадлежность кривой
      printf("Тест 8:\n");
      printf("сгенерирую случайное k в диапозоне 0 <= k <q\n");
      BIGNUM* k = BN_new();
      BN_rand_range(k, prm.q); //0<=k<q
      printf("k:%s\n", BN_bn2dec(k));
      printf("посчитаю [k]P:\n");
      CalculateDegree(&resPoint, &P_base, &curve,k);
      PrintInAffin(&resPoint, &curve);
      if(CheckPointIsOnCurve(&resPoint, &curve)){
          printf("точка [k]P находится на кривой\n\n");
      }
      else{
          printf("точка [k]P не находится на кривой\n\n");
      }


      // проверить [k1]P + [k2]P = [k1 + k2]P
      printf("Тест 9:\n");
      printf("[k1]P + [k2]P = [k1 + k2]P\n");
      BIGNUM* k1 = BN_new();
      BIGNUM* k2 = BN_new();
      printf("сгенерирую k1 и k2\n");
      BIGNUM* maxrand = BN_new();
      BN_dec2bn(&maxrand, "100000000000000000");
      BN_rand_range(k1, maxrand);
      BN_rand_range(k2, maxrand);
      printf("k1 = %s\n", BN_bn2dec(k1));
      printf("k2 = %s\n", BN_bn2dec(k2));
      BN_add(k,k1,k2);                          //k=k1+k2
      printf("k = k1 + k2 = %s\n", BN_bn2dec(k));
      struct Point res1={NULL, NULL, NULL};
      struct Point res2={NULL, NULL, NULL};
      struct Point res3={NULL, NULL, NULL};
      InitPoint(&res1,"0", "1", "1");
      InitPoint(&res2,"0", "1", "1");
      InitPoint(&res3,"0", "1", "1");
      CalculateDegree(&res1, &P_base, &curve, k1);
      CalculateDegree(&res2, &P_base, &curve, k2);
      CalculateDegree(&res3, &P_base, &curve,k);
      AddPoints(&res1, &res2, &resPoint, &curve);

      if(!(IsEqual(&resPoint, &res3, &curve))){
                printf("[k1]P + [k2]P равно [k1 + k2]P\n\n");
       }
      else{
          printf("[k1]P + [k2]P не равно [k1 + k2]P\n\n");
      }
      if(CheckPointIsOnCurve(&res3, &curve)){
          printf("точка [k]P находится на кривой\n\n");
      }
      else{
          printf("точка [k]P не находится на кривой\n\n");
      }

    FreePoint(&res3);
    FreePoint(&res2);
    FreePoint(&res1);
    BN_free(k2);
    BN_free(k1);
    BN_free(k);
    FreePoint(&negP);
    BN_free(degree);
    FreePoint(&resPoint);
    FreePoint(&P2);
    FreePoint(&P1);
    FreePoint(&E);
    FreePoint(&P_base);
    FreeParam(&prm);
    FreeJacobiCurve(&curve);

}

int main()
{
    Testing();
    return 0;
}
