#include "Utility.h"
class Myturple{
public:
    int u;
    int x;
    int y;
    int fold_vertex;
    int state;
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