// #include <bits/std++.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <thread>
#include <random>
#include <limits>
#include <string>
#include <fstream>
using  namespace std;
#define rep(i,n) for (int i = 0; i < (int)(n); i++)
using ll = long long;
using P = pair<int,int>;

// 整数を奇数部分と偶数部分に分ける、ミラーロビンに使う
vector<ll> sepEvenOld(ll N){
    vector<ll> sd={0,N};
    while(N%2==0){
        sd[0]++;
        sd[1]/=2;
        N/=2;
        // cout << sd[0] <<" " << sd[1]<< " " << N <<endl;
    }
    return sd;
};

// a^n % mod の高速計算
template<class T> T pow_mod(T a, T n, T mod){
    T res = 1;
    a %= mod;
    while (n){
        if (n&1){
            res = (res*a) % mod;
        }
        a=(a*a)%mod;
        n = n >> 1;
    }
    return res;

}

// ミラーロビン
bool MillerRabin(ll p, vector<ll> A){
    vector<ll> sd;
    sd = sepEvenOld(p-1);
    ll s=sd[0];
    ll d=sd[1];
    bool primeFlag = true;
    // cout << "p:"<< p << "  p-1=2^" << s << "*" << d <<endl;

    if (p==2){
        // cout << p << " is prime" << endl;
        return true;
    }
    if (p%2==0 || p==1){
        // cout << p << " is not prime" << endl;
        return false;
    }
    // cout <<" p="<<p<<" s="<<s<<" d="<<d<<endl; 

    for (auto a:A){
        if (p<=a) continue;
        ll t,x=pow_mod<__int128_t>(a,d,p);
        if (x!=1){
            for (t=0;t<s;++t){
                if (x==p-1)break;
                x=__int128_t(x)*x%p;
            }
            if (t==s){
                // cout << p << " is not prime" << endl;
                primeFlag = false;
                break;
            }
        }
    }
    // cout << p << " is  prime" << endl;
    return primeFlag;
}
// 素数判定
bool is_prime(ll p){
    if (p<4759123141LL){
        return MillerRabin(p,{2,7,61});
    }else{
        return MillerRabin(p,{2,325,9375,28178,450775,9780504,1795265022});
    }
}
// 最大公約数を計算
int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

ll find_prime_factor(ll n) {
    if (n % 2 == 0)  // nが偶数なら2を返す
        return 2;

    ll m = static_cast<ll>(std::pow(n, 0.125)) + 1;  // O(n^(1/8))

    for (ll c = 1; c < n; ++c) {
        auto f = [n, c](ll a) {
            return (a * a % n + c) % n;  // 疑似乱数生成関数
        };
        
        ll y = 0, g = 1, q = 1, r = 1, k = 0;
        ll x = 0, ys = 0;  // xとysをスコープの外に定義

        while (g == 1) {
            x = y;  // xをスコープ外で定義しているので問題なし
            while (k < 3 * r / 4) {  // k < 3r/4の間、GCD計算をスキップ
                y = f(y);
                ++k;
            }
            
            while (k < r && g == 1) {
                ys = y;  // バックトラック用にysを保存
                for (ll i = 0; i < std::min(m, r - k); ++i) {
                    y = f(y);
                    q = q * std::abs(x - y) % n;
                }
                g = gcd(q, n);
                k += m;
            }
            k = r;
            r *= 2;  // GCD > 1が見つからなければ、rを2倍にして繰り返す
        }

        if (g == n) {  // GCDがnの場合、バックトラックして1つずつ検証
            g = 1;
            y = ys;
            while (g == 1) {
                y = f(y);
                g = gcd(std::abs(x - y), n);
            }
        }

        if (g == n)  // GCDがnであれば、cを変えて再試行
            continue;
        
        if (is_prime(g))  // GCDが素数であればそれを返す
            return g;
        else if (is_prime(n / g))  // n / GCDが素数であればそれを返す
            return n / g;
        else  // 再帰的に素因数探索を続ける
            return find_prime_factor(g);
    }
    return n;  // 素因数が見つからない場合、nを返す
}

ll findPrimitiveRoot(ll p, vector<ll> factor){
    ll g = 1;
    bool flag = true;
    while(flag){
        g++;
        flag = false;
        for (auto f:factor){
            if (pow_mod<ll>(g,(p-1)/f,p)==1){
                flag = true;
                break;
            }
        }
    }
    return g;
}

void valPrint(ll y, ll g, ll x, ll p){     
    cout << y << " = " << g << "^" << x << " mod " << p << endl;
}

#define PRIME_RANGE 1000000

ll primeGen(){
    mt19937_64 engine(random_device{}());
    uniform_int_distribution<ll> dist(
        2,
        // numeric_limits<ll>::max()
        PRIME_RANGE
        
    );
    ll p = dist(engine);
    while(!is_prime(p)){
        p = dist(engine);
    }
    return p;
}
ll primitiveRootGen(ll p){
    ll g;
    vector<ll> factor;
    ll n=p-1;

    while (n>1){
        ll f = find_prime_factor(n);
        // cout << "f:" << f << endl;
        factor.push_back(f);
        n = n/f;
        
    }
    cout << "p-1 = ";
    for (auto f:factor){
        cout << f << ' ';
    }
    cout << endl;
    // 生成元（最小原子根）を探す
    g = findPrimitiveRoot(p,factor);
    rep(i,p-1){ 
        // valPrint(pow_mod<ll>(g,i+1,p),g,i+1,p); 
    }
    cout << "PrimitiveRoot:" << g << endl;
    return g;
}
ll bruteForce(ll y, ll g, ll p){
    ll x = 1;
    while(1){
        if (y==pow_mod<ll>(g,x,p)){ return x; }
        x++;
    }
}
ll BSGS(ll y, ll g, ll p){
    ll x=0;
    double m = ceil(sqrt(p));
    vector<ll> gt(p);
    vector<ll> gt_index;
    ll gr = 1;
    rep(i,m){
        gt[gr] = i;
        gt_index.push_back(gr);
        // cout << "gr" << gr << endl;
        gr = (gr*g)%p;
    }

    ll gm = pow_mod<ll>(g,p-m-1,p);

    ll ygqm = y; //ygqm=y*g^(-qm)
    rep(q,m){
        if (find(gt_index.begin(), gt_index.end(),ygqm) != gt_index.end()){
            return q*m + gt[ygqm];
        }
        ygqm = (ygqm*gm) % p;
    }
    return 0;
}
ll f_function(ll x, ll y, ll g, ll p){
    if (x%3==0){return (y*x)%p;}
    if (x%3==1){return (x*x)%p;}
    if (x%3==2){return (g*x)%p;}
}
ll g_function(ll x, ll a, ll p){
    if (x%3==0){return a    %((p-1)/2);}
    if (x%3==1){return (2*a)%((p-1)/2);}
    if (x%3==2){return (a+1)%((p-1)/2);}
}
ll h_function(ll x, ll b, ll p){
    if (x%3==0){return (b+1)%((p-1)/2);}
    if (x%3==1){return (2*b)%((p-1)/2);}
    if (x%3==2){return b    %((p-1)/2);}
}
ll mod_inverse(ll a, ll m) {
    ll m0 = m, t, q;
    ll x0 = 0, x1 = 1;
    if (m == 1){
        return 0;
    }
    while (a > 1) {
        cout << a << " " << m << endl;
        q = a / m;
        t = m;
        m = a % m, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0){
        x1 += m0;
    }
    return x1;
}
ll mod_inv(ll a, ll mod){
    ll inv = 1;
    rep(i,mod){
        inv++;
        // cout << "a*inv % mod: " << (a*inv) % mod << endl;
        if((a*inv) % mod ==1){ return inv; }
    }
    return -1;
}

// 拡張ユークリッドの互除法
ll gcdExtended(ll a, ll b, ll &x, ll &y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    ll x1, y1;
    ll gcd = gcdExtended(b % a, a, x1, y1);
    
    x = y1 - (b / a) * x1;
    y = x1;
    
    return gcd;
}
// modの逆元を求める関数（pは素数でなくてもOK）
ll modInverse(ll a, ll p) {
    ll x, y;
    ll gcd = gcdExtended(a, p, x, y);
    
    // 逆元が存在する場合のみ（aとpが互いに素な場合）
    if (gcd != 1) {
        // std::cerr << "Inverse doesn't exist (a and p are not coprime)" << std::endl;
        return -1;
    } else {
        // xは負の値になることがあるので、正のmodの範囲に修正
        return (x % p + p) % p;
    }
}
ll modPow(ll a, ll b ,ll mod){
    ll result = 1;
    a %= mod; // aをmodで小さくする
    while (b > 0) {
        if (b % 2 == 1) { // bが奇数なら
            result = (result * a) % mod;
        }
        a = (a * a) % mod; // aを2乗してmodを取る
        b /= 2; // bを半分にする
        // cout << result << endl;
    }
    return result;
}

// ポラード・ロー法による離散対数の解法
ll pollardRho(ll y, ll g, ll p) {
    // 初期化
    mt19937_64 engine(random_device{}());
    uniform_int_distribution<ll> dist(1,10);

    // ll x = dist(engine), a = dist(engine), b = dist(engine);
    // cout << x << a << b << endl;
    ll x = 1, a = 0, b = 0;
    ll X = x, A = a, B = b;
    ll q = (p-1)/2;

    for (int i = 1; i < p; i++) {
        // f(x)を用いてxi, ai, biの更新
        if (x % 3 == 0) {
            // x = (x * x) % p;
            x = (x * x) % p;
            a = (a * 2) % q;
            b = (b * 2) % q;
        } else if (x % 3 == 1) {
            x = (x * g) % p;
            // x = (x * g);
            a = (a + 1) % q;
        } else {
            x = (x * y) % p;
            // x = (x * y);
            b = (b + 1) % q;
        }

        // f(X)を用いてX, A, Bの更新（倍速）
        for (int j = 0; j < 2; j++) {
            if (X % 3 == 0) {
                X = (X * X) % p;
                A = (A * 2) % q;
                B = (B * 2) % q;
            } else if (X % 3 == 1) {
                X = (X * g) % p;
                A = (A + 1) % q;
            } else {
                X = (X * y) % p;
                B = (B + 1) % q;
            }
        }
        
        
        // xとXが一致した場合
        if (x == X) {
            // cout << "x:" << x << "  a:" << a << "  b:" << b << endl;
            // cout << "X:" << X << "  A:" << A << "  B:" << B << endl;
            ll r = (a - A) % q;
            // cout << "r:" << r << endl;
            if (r == 0) {
                // cout << "failed (a-A=0)" << endl;
                continue;
                // return -1; // 解けない時
            }
            while(r<0){ r += q; } //正になるまで加算
            // cout << "r:" << r << endl;
            
            ll inv = modInverse((B - b) % q, q);
            if (inv == -1){ continue;}
            // cout << "inv:" << inv << endl;
            while(inv<0){ inv += q; }
            // cout << "inv:" << inv << endl;

            if (inv == 0) {
                // cout << "not exist modular inverse " << endl;
                continue;
                // return -1;
            }
            // if (inv == -1 ){
            //     cout << "not inv" << endl;
            //     return -1;
            // }
            ll res = (r * inv) % q;
            if (modPow(g,res,p) == y){ return res;}
            if (modPow(g,res+q,p) == y){ return res + q;}
        }
    }
    cout << "not find" << endl;
    return -1;
}

int main(){
    // vector<ll> P={101,109,100022,343454543,104729,1099511627791,223372036854775807};
    int num = 3;
    vector<ll> P(num);
    vector<ll> Y(num);
    vector<ll> G(num);
    cout << "test output" << endl;
    // ll tmp = primitiveRootGen(947);

    
    rep(i,100){
        // y,g（生成元）,p（素数）の生成
        ll p = primeGen();
        // p = 19;
        P.push_back(p);

        mt19937_64 engine(random_device{}());
        uniform_int_distribution<ll> dist(2,p-1);
        ll y = dist(engine);
        // y = 7;
        Y.push_back(y);
        cout << "p=" << p << " y=" <<  y << endl;
        ll g = primitiveRootGen(p);
        G.push_back(g);

        // // ロー法
        cout << "\nPollarRho" << endl;
        auto start = chrono::high_resolution_clock::now();
        auto ans = pollardRho(y,g,p);
        auto end = chrono::high_resolution_clock::now();
        if (ans==-1){
            std::ofstream file("PollarRho.csv", std::ios::app);
            file << p << "," << -1 << endl;    
        }else{
            auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
            cout << "Exe Time is : " << duration.count() << " micro seconds " << endl;
            std::ofstream filePollar("PollarRho.csv", std::ios::app);
            filePollar << p << "," << duration.count() << endl;

            duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            std::string quotation = std::to_string(duration.count());
            cout << "Exe Time is : " << duration.count()/1000 << " ms " << endl;
            valPrint(y,g,ans,p);
            
            
            
        }
        



        // 総当たり法
        cout << "\nbruteForce" << endl;
        ll cnt = 0;
        start = chrono::high_resolution_clock::now(); // 開始時間を取得
        ans = bruteForce(y,g,p);
        end = chrono::high_resolution_clock::now();// 終了時間を取得

        auto dur = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "Exe Time is : " << dur.count() << " micro seconds " << endl;
        std::ofstream filebrute("bruteForce.csv", std::ios::app);
        std::string quo = std::to_string(dur.count()/1000);
        filebrute << p << "," << dur.count() << endl;

        dur = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Exe Time is : " << dur.count()/1000 << " ms " << endl;
        valPrint(y,g,ans,p);
        
        // // BSGS法
        cout << "\nBSGS" << endl;
        start = chrono::high_resolution_clock::now();
        ans = BSGS(y,g,p);
        end = chrono::high_resolution_clock::now();

        dur = chrono::duration_cast<chrono::microseconds>(end - start);
        cout << "Exe Time is : " << dur.count() << " micro seconds " << endl;
        std::ofstream fileBSGS("BSGS.csv", std::ios::app);
        quo = std::to_string(dur.count()/1000);
        fileBSGS << p << "," << dur.count() << endl;

        dur = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Exe Time is : " << dur.count()/1000 << " ms " << endl;
        valPrint(y,g,ans,p);
        
    }
    
    return 0;
}