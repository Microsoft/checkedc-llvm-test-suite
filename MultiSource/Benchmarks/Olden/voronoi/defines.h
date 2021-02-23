/* For copyright information, see olden_v1.0/COPYRIGHT */

#include <math_checked.h>

#pragma BOUNDS_CHECKED ON

typedef int BOOLEAN;
typedef unsigned long      uptrint;

struct edge_rec {
  _Ptr<struct VERTEX> v;
  _Array_ptr<struct edge_rec> next : quad_bounds(next);
  long wasseen;
  _Array_ptr<void> Buffer;
};


struct get_point {
  _Ptr<struct VERTEX> v;
  double curmax;
  int seed;
};


typedef _Ptr<struct edge_rec> EDGE_PTR;
typedef _Ptr<struct VERTEX> VERTEX_PTR;
typedef _Array_ptr<struct edge_rec> QUAD_EDGE;
  
struct VEC2 {
  double x,y;
  double norm;
};


struct VERTEX {
  struct VEC2 v;
  _Ptr<struct VERTEX> left, right;
};


typedef struct {
  QUAD_EDGE left : quad_bounds(left), right : quad_bounds(right);
} EDGE_PAIR;


struct EDGE_STACK {
    int ptr;
    _Array_ptr<QUAD_EDGE> elts : count(stack_size * 2);
    int stack_size;
};


#define onext(a) (a)->next
#define oprev(a) rot(onext(rot(a)))
#define lnext(a) rot(onext(rotinv(a)))
#define lprev(a) sym(onext(a))
#define rnext(a) rotinv(onext(rot(a)))
#define rprev(a) onext(sym(a))
#define dnext(a) sym(onext(sym(a)))
#define dprev(a) rotinv(onext(rotinv(a)))

#define X(r) r->v.x
#define Y(r) r->v.y
#define NORM(r) r->v.norm
#define orig(a) (a)->v
#define dest(a) orig(sym(a))
#define seen(a) (a)->wasseen
  
#define origv(a) orig(a)->v
#define destv(a) dest(a)->v


#define SIZE ((uptrint) sizeof(struct edge_rec))
#define ANDF ((uptrint) 4*sizeof(struct edge_rec) - 1)

#define sym(a) ((QUAD_EDGE) (((uptrint) (a)) ^ 2*SIZE))
#define rot(a) ((QUAD_EDGE) ( (((uptrint) (a) + 1*SIZE) & ANDF) | ((uptrint) (a) & ~ANDF) ))
#define rotinv(a) ((QUAD_EDGE) ( (((uptrint) (a) + 3*SIZE) & ANDF) | ((uptrint) (a) & ~ANDF) ))
#define base(a) ((QUAD_EDGE) (((uptrint) (a)) & ~ANDF))

#define quad_bounds(a) bounds(base(a), base(a) + 4)

QUAD_EDGE alloc_edge(void) : count(4);
void free_edge(QUAD_EDGE e : quad_bounds(e));
QUAD_EDGE makeedge(VERTEX_PTR origin, VERTEX_PTR destination) : count(4);
void splice(QUAD_EDGE a : quad_bounds(a), QUAD_EDGE b : quad_bounds(b));
void swapedge(QUAD_EDGE e : quad_bounds(e));
void deleteedge(QUAD_EDGE e : quad_bounds(e));
QUAD_EDGE build_delaunay_triangulation(VERTEX_PTR tree, VERTEX_PTR extra) : count(4);
EDGE_PAIR build_delaunay(VERTEX_PTR tree, VERTEX_PTR extra);
EDGE_PAIR do_merge(QUAD_EDGE ldo : quad_bounds(ldo), QUAD_EDGE ldi : quad_bounds(ldi), QUAD_EDGE rdi : quad_bounds(rdi), QUAD_EDGE rdo : quad_bounds(rdo));
QUAD_EDGE connect_left(QUAD_EDGE a : quad_bounds(a), QUAD_EDGE b : quad_bounds(b)) : count(4);
QUAD_EDGE connect_right(QUAD_EDGE a : quad_bounds(a), QUAD_EDGE b : quad_bounds(b)) : count(4);

int myrandom(int seed);
void zero_seen(_Ptr<struct EDGE_STACK> my_stack, QUAD_EDGE edge : quad_bounds(edge));
QUAD_EDGE pop_edge(_Ptr<struct EDGE_STACK> x) : count(4);

#define drand(seed) (((double) (seed=myrandom(seed))) / (double) 2147483647)

extern _Ptr<VERTEX_PTR> vp ;
extern _Array_ptr<struct VERTEX> va ;
extern _Ptr<EDGE_PTR> next ;
extern _Ptr<VERTEX_PTR> org ;
extern int num_vertices, num_edgeparts, stack_size ;
extern int to_lincoln , to_off, to_3d_out, to_color , voronoi , delaunay , interactive , ahost ;
extern _Array_ptr<char> see;
extern int NumNodes, NDim;


void push_ring(_Ptr<struct EDGE_STACK> stack, QUAD_EDGE edge : quad_bounds(edge));
void push_edge(_Ptr<struct EDGE_STACK> stack, QUAD_EDGE edge : quad_bounds(edge));
BOOLEAN ccw(VERTEX_PTR a, VERTEX_PTR b, VERTEX_PTR c);
_Unchecked
int dealwithargs(int argc, _Array_ptr<char *> argv : count(argc));
void output_voronoi_diagram(QUAD_EDGE edge : quad_bounds(edge), int nv, _Ptr<struct EDGE_STACK> stack);
struct get_point get_points(int n, double curmax,int i, int seed,
                            int processor, int numnodes);
_Ptr<struct EDGE_STACK> allocate_stack(int num_vertices);


double V2_cprod(struct VEC2 u,struct VEC2 v);
double V2_dot(struct VEC2 u, struct VEC2 v);
struct VEC2 V2_times(double c, struct VEC2 v);
struct VEC2 V2_sum(struct VEC2 u, struct VEC2 v);
struct VEC2 V2_sub(struct VEC2 u, struct VEC2 v);
double V2_magn(struct VEC2 u);
struct VEC2 V2_cross(struct VEC2 v);

#pragma BOUNDS_CHECKED OFF
