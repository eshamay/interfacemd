#include <iostream>
#include <vector>

using namespace std;


struct binner : public unary_function<int,int>
{
    typedef std::vector<int> binner_t;
    binner_t values;

    void operator() (const int i)
    {   
        values.push_back(i);
    }   
    
    binner_t& Values () { return values; }
    int size() { return int(values.size()); }
};

int main() {

    binner bins;
    for (int i = 0; i < 20; i++)
        bins(int(rand()*20.0/RAND_MAX));

    for (int i = 0; i < bins.size(); i++)
        printf ("%d\n", bins.Values()[i]);

  return 0;
}

