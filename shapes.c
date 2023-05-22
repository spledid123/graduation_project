// 求各种形状的SDF,闵可夫斯基和

#include "svpng.inc"
#include <math.h>   // fminf(), sinf(), cosf(), sqrt()
#include <stdlib.h> // rand(), RAND_MAX

#define TWO_PI 6.28318530718f
#define W 512
#define H 512
#define N 64
#define MAX_STEP 64
#define MAX_DISTANCE 2.0f
#define EPSILON 1e-6f
#define DELTA 0.1f

typedef struct
{
    double sd, emissive;
} Result;
// 定义物体

unsigned char img[W * H * 3];

double circleSDF(double x, double y, double cx, double cy, double r)
{
    double ux = x - cx, uy = y - cy;
    return sqrtf(ux * ux + uy * uy) - r;
}

double planeSDF(double x, double y, double px, double py, double nx, double ny)
{
    return (x - px) * nx + (y - py) * ny;
}
// 平面SDF,p为平面内一点，n为单位矢量，x=(x,y)为任意一点
// SDF=(x-p)·n

double segmentSDF(double x, double y, double ax, double ay, double bx, double by)
{
    double vx = x - ax, vy = y - ay, ux = bx - ax, uy = by - ay;
    double t = fmaxf(fminf((vx * ux + vy * uy) / (ux * ux + uy * uy), 1.0f), 0.0f);
    double dx = vx - ux * t, dy = vy - uy * t;
    return sqrtf(dx * dx + dy * dy);
}
// 直线
double segmentSDF2(double x, double y, double ax, double ay, double bx, double by)
{
    double cross = (bx - ax) * (x - ax) + (by - ay) * (y - ay);
    if (cross <= 0)
        return sqrtf((x - ax) * (x - ax) + (y - ay) * (y - ay));

    double d2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
    if (cross >= d2)
        return sqrtf((x - bx) * (x - bx) + (y - by) * (y - by));

    double r = cross / d2;
    double px = ax + (bx - ax) * r;
    double py = ay + (by - ay) * r;
    return sqrtf((x - px) * (x - px) + (py - y) * (py - y));
}
//直线，据说是穷人版

double capsuleSDF(double x, double y, double ax, double ay, double bx, double by, double r)
{
    return segmentSDF(x, y, ax, ay, bx, by) - r;
}
//胶囊形，a,b分别为两圆心坐标，r为半径
//SDF=点x到ab所在直线u的距离+r

double boxSDF(double x, double y, double cx, double cy, double theta, double sx, double sy)
{
    double costheta = cosf(theta), sintheta = sinf(theta);
    double dx = fabs((x - cx) * costheta + (y - cy) * sintheta) - sx;
    double dy = fabs((y - cy) * costheta - (x - cx) * sintheta) - sy;
    double ax = fmaxf(dx, 0.0f), ay = fmaxf(dy, 0.0f);
    return fminf(fmaxf(dx, dy), 0.0f) + sqrtf(ax * ax + ay * ay);
}
//矩形SDF，中心点c，旋转角\theta(方块固连坐标系与空间坐标系的角度),
//半对角线长s
//先算出矩形体坐标系x的坐标x'，x=cot(\theta)x'+c,则x'=cot(-\theta)(x-c)
//x'=(x - cx) * costheta + (y - cy) * sintheta,y'=(y - cy) * costheta - (x - cx) * sintheta
// dx=|x'|-sx;dy=|y'|-sy;
//计算方法很妙，自行体会dx,dy分别大于小于等于0的情况，SDF是怎么算的

double triangleSDF(double x, double y, double ax, double ay, double bx, double by, double cx, double cy)
{
    double d = fminf(fminf(
                        segmentSDF(x, y, ax, ay, bx, by),
                        segmentSDF(x, y, bx, by, cx, cy)),
                    segmentSDF(x, y, cx, cy, ax, ay));
    return (bx - ax) * (y - ay) > (by - ay) * (x - ax) &&
                   (cx - bx) * (y - by) > (cy - by) * (x - bx) &&
                   (ax - cx) * (y - cy) > (ay - cy) * (x - cx)
               ? -d
               : d;
}
double clamp(double x, double xmin, double xmax)
{
    if (xmin > xmax)
    {
        printf("clamp is wrong!");
        return 0;
    }
    if (x < xmin)
        return xmin;
    if (x > xmax)
        return xmax;
    return x;
}
int sign(double x)
{
    if (x > 0.0f)
        return 1;
    if (x < 0.0f)
        return -1;
    return 0;
}
//double sdHexagon( in vec2 p, in double r )
double num_atanf(double px, double py)
{
    double PI = 3.14159265f;
    if (px == 0)
    {
        if (py == 0)
            return 0;
        if (py > 0)
            return PI / 2;
        if (py < 0)
            return -PI / 2;
    }
    if (px > 0)
    {
        return atan(py / px);
    }
    if (px < 0)
    {
        return atan(py / px) + PI;
    }
}
double myNPolygonSDF(double x, double y, double cx, double cy, double r, int n)
//double myNPolygonSDF( in vec2 p, in double r, in double N )
{
    //预计算一些计算需要的值
    double an = 6.2831853 / (double)n;
    double he = r * tan(0.5 * an) - DELTA;
    double px = x - cx;
    double py = y - cy;

    double bn;
    double rx = r - DELTA, ry = (r - DELTA) * tan(0.5 * an); //小圆圆心
    //这句代码是将图形关于y=-x直线做对称，可以使正多边形顶角向上
    //p = -p.yx;
    bn = -an * floor((num_atanf(px, py) + 0.5 * an) / an);
    //printf("bn=%f\n", bn / 3.14159 * 180);
    //旋转至①区域
    //printf("px=%f,py=%f,atan=%f,px*px+py*py=%f\n", px, py, (num_atanf(px, py)) / 3.14159 * 180, px * px + py * py);
    double csx = cos(bn);
    double csy = sin(bn);
    //printf("csx=%f,csy=%f\n", csx, csy);
    //vec2  cs = vec2(cos(bn),sin(bn));
    //这里注意glsl里矩阵列优先储存
    double pxx=px;
    px = (px * csx) - (py * csy);
    py = py * csx + pxx * csy;
    //printf("px=%f,py=%f,px*px+py*py=%f,cla=%f,he=%f,atan=%f,an/2=%f\n", px, py, px * px + py * py, clamp(py, -he, he), he, num_atanf(px, py) / 3.14159 * 180, 0.5 * an / 3.14159 * 180);
    //p = mat2(cs.x,cs.y,-cs.y,cs.x)*p;
    //①区域的SDF
    if (0 < num_atanf(px - rx, py - ry))
    {
        if (num_atanf(px - rx, py - ry) < 0.5 * an)
            return sqrt((px - rx) * (px - rx) + (py - ry) * (py - ry)) - DELTA;
    }
    if (0 > num_atanf(px - rx, py + ry))
    {
        if (num_atanf(px - rx, py + ry) + 0.5 * an > 0)
            return sqrt((px - rx) * (px - rx) + (py + ry) * (py + ry)) - DELTA;
    }

    px = px - r;
    py = py - clamp(py, -he, he);

    double lp = sqrt(px * px + py * py);
    return lp * sign(px);
    //return length(p-vec2(r,clamp(p.y,-he,he)))*sign(p.x-r);
}

Result intersectOp(Result a, Result b)
{
    return a.sd > b.sd ? a : b;
}
// 交集

Result scene(double x, double y)
{
    /*Result a = {circleSDF(x, y, 0.5f, 0.5f, 0.2f), 1.0f};
    Result b = {planeSDF(x, y, 0.0f, 0.5f, 0.0f, 1.0f), 0.8f};
    Result c = {capsuleSDF(x, y, 0.4f, 0.4f, 0.6f, 0.6f, 0.1f), 1.0f};*/
    Result d = {boxSDF(x, y, 0.5f, 0.5f, TWO_PI / 16.0f, 0.3f, 0.1f), 1.0f};
    /*
    Result e = {boxSDF(x, y, 0.5f, 0.5f, TWO_PI / 16.0f, 0.3f, 0.1f) - 0.1f, 1.0f};
    Result f = {triangleSDF(x, y, 0.5f, 0.2f, 0.8f, 0.8f, 0.3f, 0.6f), 1.0f};
    Result g = {triangleSDF(x, y, 0.5f, 0.2f, 0.8f, 0.8f, 0.3f, 0.6f) - 0.1f, 1.0f};
   */
    //Result h = {myNPolygonSDF(x, y, 0.5f, 0.5f, 0.2f, 4), 1.0f};

    // return a;
    // return b;
    // return intersectOp(a, b);
    // return c;
    // return d;
    // return e;
    // return f;
    // return g;
    return d;
}

double trace(double ox, double oy, double dx, double dy)
{
    double t = 0.0f;
    for (int i = 0; i < MAX_STEP && t < MAX_DISTANCE; i++)
    {
        Result r = scene(ox + dx * t, oy + dy * t);
        if (r.sd < EPSILON)
            return r.emissive;
        t += r.sd;
    }
    return 0.0f;
}

double sample(double x, double y)
{
    double sum = 0.0f;
    for (int i = 0; i < N; i++)
    {
        double a = TWO_PI * (i + (double)rand() / RAND_MAX) / N;
        sum += trace(x, y, cosf(a), sinf(a));
    }
    return sum / N;
}

int main()
{
    unsigned char *p = img;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++, p += 3)
            p[0] = p[1] = p[2] = (int)(fminf(sample((double)x / W, 1.0f - (double)y / H) * 255.0f, 255.0f));
    svpng(fopen("shapes(2).png", "wb"), W, H, img, 0);
}