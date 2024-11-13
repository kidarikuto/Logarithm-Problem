def pollard_rho(g, y, p):
    q = (p-1) // 2
    # ランダムに移動するための関数
    def new_xab(x, a, b,  g, y, p, q):
        subset = x % 3
        if subset == 0:
            return ((x*x) % p, (a*2) % q, (b*2) % q)
        if subset == 1:
            return ((x*g) % p, (a+1) % q, b        )
        if subset == 2:
            return ((x*y) % p, a        , (b+1) % q)
    # フロイドの循環検出法
    x, a, b = 1, 0, 0
    X, A, B = x, a, b
    for i in range(1, p):
        x, a, b = new_xab(x, a, b,  g, y, p, q)
        X, A, B = new_xab(X, A, B,  g, y, p, q)
        X, A, B = new_xab(X, A, B,  g, y, p, q)
        print(f"{x=} {a=} {b=} {X=} {A=} {B=}")
        if x == X:
            break
    res = ((a - A) * pow(B - b, -1, q)) % q
    print(pow(B - b, -1, q))
    print(f"{res=}")

    if pow(g, res, p) == y:
        return res
    if pow(g, res + q, p) == y:
        return res + q
    return None

g = 2
y = 7
p = 19
x = pollard_rho(g, y, p)
print(x)
print(pow(g, x, p) == y)