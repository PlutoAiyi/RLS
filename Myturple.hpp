#include "Utility.h"
class Myturple{
public:
    int u;
    int x;
    int y;
    int fold_vertex;
    int state;
    //state=-1表示处理过了
    //state=0代表u和x互斥，两个点有且只有一个在解集中
    //具体操作：从S中检测x在不在来判断是否加u
    //state=1代表xy都不在时，u一定在，x和y某个在时，u不在(x和y互斥)
    //具体操作：从S中检测x和y是否存在来判断是否加u（c,h,j,m）
    //state=2代表x和y有一个不在时，u在，若xy都在，则u不在
    //具体操作：从S中检测x和y是否都存在来判断是否加u（d）


    //// state=3代表若fold的顶点v'被纳入了，则xy都在，否则u在
    //// 具体操作：从S中检测v'在不在来判断是否加u或者xy


    //state=3现在代表x不在或y在时，u就在，只有x在且y不在的情况下u不在



    //state=4代表fold的顶点v'被纳入了，则ux都在，否则都不在


    //state=5代表y在且x不在则u在，y不在或者x在则u不在

    //state=6代表x在或者y在则u在

    //state=7代表x在且y在则u在
    Myturple(ui u,ui x,ui y,ui fold_vertex,int state){
        this->u=u;
        this->x=x;
        this->y=y;
        this->fold_vertex=fold_vertex;
        this->state=state;
    }
    Myturple(ui u,ui x,ui y,int state){
        this->u=u;
        this->x=x;
        this->y=y;
        this->state=state;
        this->fold_vertex=-1;
    }
    Myturple(ui u,ui x,int state){
        this->y=-1;
        this->fold_vertex=-1;
        this->u=u;
        this->x=x;
        this->state=state;
    }
    Myturple(Myturple* myturple){
        this->u=myturple->u;
        this->x=myturple->x;
        this->y=myturple->y;
        this->fold_vertex=myturple->fold_vertex;
        this->state=myturple->state;
    }
    Myturple(){

    }
};