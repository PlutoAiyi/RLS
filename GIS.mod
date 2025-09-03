/*********************************************
 * OPL 22.1.1.0 Model
 * Author: timelapse
 * Creation Date: 27/08/2025 at 6:37:48 PM
 *********************************************/
execute {
	cplex.epgap = 0.0;   
    cplex.tilim = 100000;      // time limit
    cplex.mipdisplay = 4;    // show intermediate solutions
	cplex.threads = 1;         // restrict solver to 1 thread
    cplex.parallelmode = -1;   // disable parallelism (optional, -1 = automatic, 0 = opportunistic, 1 = deterministic)
}








tuple edge {
    int u;
    int v;
}

 // nodes
int     n       = ...;
range   nodes  = 1..n;
int     weight[nodes] = ...;

// Edges -- sparse set
//tuple       edge        {int i; int j;}
{edge} 		Edges      =...;
int         penalty[Edges] = ...;
setof(edge) HardEdges       =...;


// Decision variables
dvar boolean x[nodes];
dvar boolean y[Edges];


/*****************************************************************************
 *
 * MODEL
 * 
 *****************************************************************************/

// Objective
maximize sum (i in nodes) weight[i]*x[i]- sum (e in Edges) penalty[e]*y[e];
subject to {
    forall(e in HardEdges)
        x[e.u] + x[e.v] <= 1;

    forall(e in Edges) {
        x[e.u] + x[e.v] - y[e] <= 1;
        x[e.u] - y[e] >= 0;
        x[e.v] - y[e] >= 0;
    }
}


// Post-processing: executed after solve
execute {
	var pyList = "[";
   	for (var i in nodes) {
   	   if (i > 0) pyList += ", ";
   	   pyList += x[i];   // <-- decision variable value (0 or 1)
   	}
   	pyList += "]";
   	writeln("x = ", pyList);
    var sumX = 0;
    for(var i in nodes) {
        sumX += x[i];
    }
    writeln("Final sum of x = ", sumX);
}


