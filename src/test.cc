#include<iostream>  
#include<cstring>  
#include<cstdlib>  
#include<cstdio>  
#include<climits>  
#include<algorithm>  
#include<fstream>
using namespace std;  

const int N = 1000;  
const int INF = 0xffffff;  
int w[N][N];//权值  
int lx[N],ly[N]; //顶标  
int linky[N];//记录与i匹配的顶点  
int visx[N],visy[N];  
int slack[N];//松弛量  
int nx,ny;//二分图两边的顶点数  

void init()  
{  
    memset(linky,-1,sizeof(linky));//记录与i匹配的顶点  
    memset(ly,0,sizeof(ly));///初始化顶标y为0  
    for(int i = 0; i < nx; i++) {
        lx[i] = -INF;
        for(int j = 0; j < ny; j++)  
        {  
            if(w[i][j] > lx[i])  
                lx[i] = w[i][j];///初始化顶标x为与顶点Xi关联的边的最大权  
        } 
    }
}  

bool find(int x)//匈牙利算法  
{  
    visx[x] = true;  
    for(int y = 0; y < ny; y++)  
    {  
        if(visy[y])  
            continue;  
        int t = lx[x] + ly[y] - w[x][y];//若t==0，则为最大权匹配；  

        if(t==0)  
        {  
            visy[y] = true;  
            if(linky[y]==-1 || find(linky[y]))  
            {  
                linky[y] = x;  
                return true;        //找到增广轨  
            }  
        }  

        else if(slack[y] > t)  
            slack[y] = t;  
    }  
    return false;                   //没有找到增广轨（说明顶点x没有对应的匹配，与完备匹配(相等子图的完备匹配)不符）  
}  

int KM()                //返回最优匹配的值  
{  
    init();  
    for(int x = 0; x < nx; x++)  
    {  
        for(int i = 0; i < ny; i++)  
            slack[i] = INF;//松弛函数初始化为无穷大  

        while(1)  
        {  
            memset(visx,0,sizeof(visx));  
            memset(visy,0,sizeof(visy));  
            if(find(x))                     //找到增广轨，退出  
                break;  
            int d = INF;  
            for(int i = 0; i < ny; i++)          //没找到，对l做调整(这会增加相等子图的边)，重新找  
            {  
                if(!visy[i] && d > slack[i])  
                    d = slack[i];  
            }  
            for(int i = 0; i < nx; i++)//修改x的顶标  
            {  
                if(visx[i])  
                    lx[i] -= d;  
            }  
            for(int i = 0; i < ny; i++)//修改y的顶标  
            {  
                if(visy[i])  
                    ly[i] += d;  
                else  
                    slack[i] -= d;//修改顶标后，不在交错树中的y顶点的slack值都要减去d；  
            }  
        }  

    }  

    int result = 0;  
    for(int i = 0; i < ny; i++)  
    {  
        if(linky[i]>-1)  
            result += w[linky[i]][i];  
    }  
    return result;  
}


int main() {
    ifstream ifs("kd.txt");
    while (ifs >> nx >> ny) {
        if (!nx || !ny)
            break;
        int a, b, c;
        while (ifs >> a >> b >> c) {
            cout << "add node from left index " << a << " to right index " << b << " with weight " << c << endl; 
            w[a][b] = c;
        }
        printf("%d\n", KM());
    }
    return 0;
}
