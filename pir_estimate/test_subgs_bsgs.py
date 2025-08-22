'''
    Test BSGS PIR

'''

import math

# copy from test_our_pir.py
def sum_power(a, r):
    sum = 0
    for i in range(r):
        sum += pow(a, i)
    return sum

def compute_sigma2_MatrixVectorMul(N, p, ell, Bg, P, sigma = 8):
    rtq = 0
    sigma2_mul = N * N * pow(p/2., 1) * (2 * math.pi) + rtq

    sigma2_rotate = (1 + N * ell * pow(Bg, 2) / 4) * pow(sigma, 2)
    sigma2_scale = 2 * math.pi * (round(2 * N / 3.) + 1) * pow(P/2, 2)

    # sigma2_acc = sqrtN * pow(math.floor(p / 2), 2) * pow(sigma, 2)
    sigma2_total = N / 2 * N * pow(p, 2) * (sigma2_rotate + sigma2_scale) / 4

    return sigma2_mul + sigma2_total

def compute_sigma2_bsgs(N, p, ell, Bg, P, ell_gs, Bg_gs, P_gs, N1, h = 2731, sigma = 8, q = pow(2, 56), Q = pow(2, 80)):
    rtq = 0
    sigma2_mul = N * N * pow(p/2., 1) * (2 * math.pi) + rtq

    # the factor 2 is for two phase for hoisting
    sigma2_rotate = (1 + 2 * N * ell * pow(Bg, 2) / 4) * pow(sigma, 2)
    sigma2_scale = 2 * math.pi * (h + 1) * pow(P/2, 2)

    sigma2_modswitch = 2 * math.pi * (h + 1) / 2
    Qq = Q / q
    sigma2_after_modswitch = (sigma2_rotate + sigma2_scale) / Qq / Qq + sigma2_modswitch

    # sigma2_acc = sqrtN * pow(math.floor(p / 2), 2) * pow(sigma, 2)
    sigma2_total = N1 * N * pow(p, 2) * (sigma2_after_modswitch) / 4

    sigma2_gs_acc = N * ell_gs * pow(Bg_gs, 2) / 4 * pow(sigma, 2)
    sigma2_gs_scale = 2 * math.pi * (h + 1) * pow(P_gs/2, 2)

    N2 = N / 2 / N1
    sigma2_gs_total = N2 * sigma2_total + (N2 - 1) * (sigma2_gs_acc + sigma2_gs_scale)

    return sigma2_mul + sigma2_gs_total

def compute_sigma2_bsgs_stateless(N, p, ell, Bg, P, ell_gs, Bg_gs, P_gs, N1, h = 2731, sigma = 8, q = pow(2, 56), Q = pow(2, 80)):
    rtq = 0
    sigma2_mul = N * N * pow(p/2., 1) * (2 * math.pi) + rtq

    sigma2_rotate = (1 + N * ell * pow(Bg, 2) / 4) * pow(sigma, 2)
    sigma2_scale = 2 * math.pi * (h + 1) * pow(P/2, 2)

    sigma2_modswitch = 2 * math.pi * (h + 1) / 2
    Qq = Q / q
    # (N1 - 1) *
    sigma2_after_modswitch = (N1 - 1) * (sigma2_rotate + sigma2_scale) / Qq / Qq + sigma2_modswitch

    # sigma2_acc = sqrtN * pow(math.floor(p / 2), 2) * pow(sigma, 2)
    sigma2_total = N1 * N * pow(p, 2) * (sigma2_after_modswitch) / 4

    sigma2_gs_acc = N * ell_gs * pow(Bg_gs, 2) / 4 * pow(sigma, 2)
    sigma2_gs_scale = 2 * math.pi * (h + 1) * pow(P_gs/2, 2)

    N2 = N / 2 / N1
    sigma2_gs_total = N2 * sigma2_total + (N2 - 1) * (sigma2_gs_acc + sigma2_gs_scale)

    return sigma2_mul + sigma2_gs_total

def compute_s2_second_dimension(N, s2_first_dimension, l_exP, Bg_exP, base_exP, h = 2731, t = pow(2, 28), s2_origin = 64):
    s2_exProduct_acc = 2 * l_exP * N * pow(Bg_exP/2, 2) * s2_origin
    s2_exProduct_scale = (h + 1) * pow(base_exP/2, 2)
    s2_second_dimension = s2_first_dimension + s2_exProduct_acc + s2_exProduct_scale

    # print("sqrt(s2_second_dim) = ", math.sqrt(s2_second_dimension))
    # print("--------log p (second dim) = ", math.log2(compute_exp_t2(t, s2_second_dimension)))

    return s2_second_dimension

def compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, P, ell_gs, Bg_gs, P_gs, N1, h = 2731, sigma = 8, q = pow(2, 56), Q = pow(2, 80), QQ = pow(2, 108)):
    ell_offline = 5
    Bg_offline = pow(2, 18)
    P_offline = pow(2, 18)
    sigma2_rotate_offline = (1 + N * ell_offline * pow(Bg_offline, 2) / 4) * pow(sigma, 2)
    sigma2_scale_offline = 2 * math.pi * (h + 1) * pow(P_offline/2, 2)
    QQQ = QQ / Q
    sigma2_offline = N1 * (sigma2_rotate_offline + sigma2_scale_offline) / QQQ / QQQ

    rtq = 0
    sigma2_mul = N * N * pow(p/2., 1) * (2 * math.pi) + rtq

    sigma2_rotate = (1 + N * ell * pow(Bg, 2) / 4) * pow(sigma, 2)
    sigma2_scale = 2 * math.pi * (h + 1) * pow(P/2, 2)

    sigma2_modswitch = 2 * math.pi * (h + 1) / 2
    Qq = Q / q
    sigma2_after_modswitch = (sigma2_rotate + sigma2_scale) / Qq / Qq + sigma2_modswitch

    # sigma2_acc = sqrtN * pow(math.floor(p / 2), 2) * pow(sigma, 2)
    sigma2_total = N1 * N * pow(p, 2) * (sigma2_after_modswitch) / 4

    sigma2_gs_acc = N * ell_gs * pow(Bg_gs, 2) / 4 * pow(sigma, 2)
    sigma2_gs_scale = 2 * math.pi * (h + 1) * pow(P_gs/2, 2)

    N2 = N / 2 / N1
    sigma2_gs_total = N2 * sigma2_total + (N2 - 1) * (sigma2_gs_acc + sigma2_gs_scale)

    return sigma2_mul + sigma2_gs_total

# copy from test_our_pir.py
# use formula
def compute_packing_lwes_formula(n, N, l, Bg, P, s2_second, h = 2731, s2_origin = 64):
    r = int(math.log2(n))
    logNn = int(math.log2(N/n))

    s2_acc = l * N * pow(Bg/2, 2) * s2_origin
    s2_scale = (h + 1) * pow(P/2, 2) * 2 * math.pi
    s2_acc_sum = sum_power(4, r) * (s2_acc + s2_scale)

    s2_second_acc = pow(4, r) * s2_second
    s2_packing = s2_second_acc + s2_acc_sum
    # s2_total = pow(2, logNn) * s2_packing + sum_power(2, logNn) * (s2_acc + s2_scale)

    #return s2_total
    return s2_packing

def compute_packing_lwes_formula_new(n, N, l, Bg, P, s2_second, h = 2731, s2_origin = 64):
    r = int(math.log2(n))
    logNn = int(math.log2(N/n))

    s2_acc = l * N * pow(Bg/2, 2) * s2_origin
    s2_scale = (h + 1) * pow(P/2, 2) * 2 * math.pi
    # s2_acc_sum = sum_power(2, r) * (s2_acc + s2_scale)
    s2_acc_sum = r * (s2_acc + s2_scale)

    s2_packing = s2_second + s2_acc_sum
    # s2_total = pow(2, logNn) * s2_packing + sum_power(2, logNn) * (s2_acc + s2_scale)

    #return s2_total
    return s2_packing

def compute_modswitch(s2_first_dimension, p, q, h = 2731):
    s2_scale = (h + 1) * 0.25
    s2_scale = s2_scale * 2 * math.pi

    s2_qQ = pow(p / q, 2) * s2_first_dimension
    s2_modswitch = s2_qQ + s2_scale

    return s2_modswitch

def compute_exp(logq, logp, sigma2):
    # x = 5
    #p = 1 - math.erf(x)
    #print(p)

    t_div_sigma = pow(2, logq-logp-1) / math.sqrt(sigma2)

    b = 2 * math.exp(-math.pi * pow(t_div_sigma, 2))

    return b

def test_probability_unsafe():
    N = 256
    p = 7681
    ell = 6
    Bg = 0x01 << 7
    PP = 0x01 << 8

    sigma2_first_dim = compute_sigma2_MatrixVectorMul(N, p, ell, Bg, PP)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(50, math.log2(p), sigma2_first_dim)))

def test_probability():
    N = 4096
    p = 40961
    ell = 0
    Bg = 0 # 0x01 << 1 # 0x01 << 7
    PP = 0 # 0x01 << 8

    sigma2_first_dim = compute_sigma2_MatrixVectorMul(N, p, ell, Bg, PP)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(50, math.log2(p), sigma2_first_dim)))

def test_probability_rns():
    N = 4096
    p = pow(2, 16) # 1228900000
    ell = 8
    Bg = 0x01 << 10 # 0x01 << 7
    PP = 0x01 << 20 # 0x01 << 8

    sigma2_first_dim = compute_sigma2_MatrixVectorMul(N, p, ell, Bg, PP)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(100, math.log2(p), sigma2_first_dim)))


# 36 MB
def test_probability_bsgs_36MB():
    N = 4096
    p = 65537
    h = 2731
    # h = 256
    # h = 128

    ell = 3
    Bg = 0x01 << 20
    PP = 0x01 << 20

    N1 = 64

    ell_gs = 2 # 3
    Bg_gs = 0x01 << 18 # 0x01 << 14
    PP_gs = 0x01 << 20 # 0x01 << 14
    # ell_gs = 49 # 3
    # Bg_gs = 0x01 << 1 # 0x01 << 12
    # PP_gs = 0x01 << 1 # 0x01 << 14

    sigma2_first_dim = compute_sigma2_bsgs(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)
    # sigma2_first_dim = compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(52, math.log2(p), sigma2_first_dim)))

    # if we add a modswitch here
    modswitch_p1bits = 50
    sigma2_first_dim = compute_modswitch(sigma2_first_dim, pow(2, modswitch_p1bits), pow(2, 56), h)
    # print("log p = ", math.log2(compute_exp(modswitch_p1bits, math.log2(p), sigma2_first_dim)))

    # external product
    l_exP = 2
    Bg_exP = pow(2, 16) # pow(2, 18)
    base_exP = pow(2, 18) # pow(2, 20)
    sigma2_second_dim = compute_s2_second_dimension(N, sigma2_first_dim, l_exP, Bg_exP, base_exP, h)


    # packing lwes
    l_atm = 4 # 3 # 4
    Bg_atm = pow(2, 13) # pow(2, 14) # pow(2, 10)
    base_atm = 0 # pow(2, 14) # pow(2, 16)
    s2_total = compute_packing_lwes_formula_new(4, N, l_atm, Bg_atm, base_atm, sigma2_second_dim)

    print("s2_total = ", s2_total)
    print("log p = ", math.log2(compute_exp(47, math.log2(p), s2_total)))

    modswitch_pbits = 26
    s2_modswitch = compute_modswitch(s2_total, pow(2, modswitch_pbits), pow(2, 50), h)
    print("log p = ", math.log2(compute_exp(modswitch_pbits, math.log2(p), s2_modswitch)))
# 256 MB
def test_probability_bsgs():
    N = 4096
    p = 65537 # 40961 # 65537
    h = 2731
    # h = 256
    # h = 128

    ell = 3
    Bg = 0x01 << 20
    PP = 0x01 << 20

    N1 = 128 # 64

    ell_gs = 2 # 3
    Bg_gs = 0x01 << 18 # 0x01 << 14
    PP_gs = 0x01 << 20 # 0x01 << 14
    # ell_gs = 49 # 3
    # Bg_gs = 0x01 << 1 # 0x01 << 12
    # PP_gs = 0x01 << 1 # 0x01 << 14

    sigma2_first_dim = compute_sigma2_bsgs(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)
    # sigma2_first_dim = compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(52, math.log2(p), sigma2_first_dim)))

    # if we add a modswitch here
    modswitch_p1bits = 50
    sigma2_first_dim = compute_modswitch(sigma2_first_dim, pow(2, modswitch_p1bits), pow(2, 56), h)
    # print("log p = ", math.log2(compute_exp(modswitch_p1bits, math.log2(p), sigma2_first_dim)))

    # external product
    l_exP = 2
    Bg_exP = pow(2, 18) # pow(2, 18)
    base_exP = pow(2, 20) # pow(2, 20)
    sigma2_second_dim = compute_s2_second_dimension(N, sigma2_first_dim, l_exP, Bg_exP, base_exP)


    # packing lwes
    l_atm = 4 # 3 # 4
    Bg_atm = pow(2, 13) # pow(2, 14) # pow(2, 10)
    base_atm = 0 # pow(2, 14) # pow(2, 16)
    l_atm = 3 # 3
    Bg_atm = pow(2, 19) # pow(2, 14)
    base_atm = 0 # pow(2, 14)

    s2_total = compute_packing_lwes_formula_new(16, N, l_atm, Bg_atm, base_atm, sigma2_second_dim)

    print("s2_total = ", s2_total)
    print("log p = ", math.log2(compute_exp(47, math.log2(p), s2_total)))

    modswitch_pbits = 26
    s2_modswitch = compute_modswitch(s2_total, pow(2, modswitch_pbits), pow(2, 50), h)
    print("log p = ", math.log2(compute_exp(modswitch_pbits, math.log2(p), s2_modswitch)))
# 1 GB
def test_probability_bsgs_1GB():
    N = 4096
    p = 65537 # 40961 # 65537
    h = 2731
    # h = 256
    # h = 128

    ell = 3
    Bg = 0x01 << 20
    PP = 0x01 << 20

    N1 = 128 # 64

    ell_gs = 2 # 3
    Bg_gs = 0x01 << 18 # 0x01 << 14
    PP_gs = 0x01 << 20 # 0x01 << 14
    # ell_gs = 49 # 3
    # Bg_gs = 0x01 << 1 # 0x01 << 12
    # PP_gs = 0x01 << 1 # 0x01 << 14

    sigma2_first_dim = compute_sigma2_bsgs(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)
    # sigma2_first_dim = compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(52, math.log2(p), sigma2_first_dim)))

    # if we add a modswitch here
    modswitch_p1bits = 50
    sigma2_first_dim = compute_modswitch(sigma2_first_dim, pow(2, modswitch_p1bits), pow(2, 56), h)
    # print("log p = ", math.log2(compute_exp(modswitch_p1bits, math.log2(p), sigma2_first_dim)))

    # external product
    l_exP = 2
    Bg_exP = pow(2, 16) # pow(2, 18)
    base_exP = pow(2, 18) # pow(2, 20)
    sigma2_second_dim = compute_s2_second_dimension(N, sigma2_first_dim, l_exP, Bg_exP, base_exP)


    # packing lwes
    l_atm = 4 # 3 # 4
    Bg_atm = pow(2, 13) # pow(2, 14) # pow(2, 10)
    base_atm = 0 # pow(2, 14) # pow(2, 16)
    s2_total = compute_packing_lwes_formula_new(16 * 4, N, l_atm, Bg_atm, base_atm, sigma2_second_dim)

    print("s2_total = ", s2_total)
    print("log p = ", math.log2(compute_exp(47, math.log2(p), s2_total)))

    modswitch_pbits = 26
    s2_modswitch = compute_modswitch(s2_total, pow(2, modswitch_pbits), pow(2, 50), h)
    print("log p = ", math.log2(compute_exp(modswitch_pbits, math.log2(p), s2_modswitch)))

# 2 GB
def test_probability_bsgs_2GB():
    N = 4096
    p = 65537 # 40961 # 65537
    h = 2731
    # h = 256
    # h = 128

    ell = 3
    Bg = 0x01 << 20
    PP = 0x01 << 20

    N1 = 128 # 64

    ell_gs = 2 # 3
    Bg_gs = 0x01 << 18 # 0x01 << 14
    PP_gs = 0x01 << 20 # 0x01 << 14
    # ell_gs = 49 # 3
    # Bg_gs = 0x01 << 1 # 0x01 << 12
    # PP_gs = 0x01 << 1 # 0x01 << 14

    sigma2_first_dim = compute_sigma2_bsgs(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)
    # sigma2_first_dim = compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(52, math.log2(p), sigma2_first_dim)))

    # if we add a modswitch here
    modswitch_p1bits = 50
    sigma2_first_dim = compute_modswitch(sigma2_first_dim, pow(2, modswitch_p1bits), pow(2, 56), h)
    # print("log p = ", math.log2(compute_exp(modswitch_p1bits, math.log2(p), sigma2_first_dim)))

    # external product
    l_exP = 2
    Bg_exP = pow(2, 16) # pow(2, 18)
    base_exP = pow(2, 18) # pow(2, 20)
    sigma2_second_dim = compute_s2_second_dimension(N, sigma2_first_dim, l_exP, Bg_exP, base_exP)


    # packing lwes
    l_atm = 4 # 3 # 4
    Bg_atm = pow(2, 13) # pow(2, 14) # pow(2, 10)
    base_atm = 0 # pow(2, 14) # pow(2, 16)
    s2_total = compute_packing_lwes_formula_new(16 * 4 * 2, N, l_atm, Bg_atm, base_atm, sigma2_second_dim)

    print("s2_total = ", s2_total)
    print("log p = ", math.log2(compute_exp(47, math.log2(p), s2_total)))

    modswitch_pbits = 26
    s2_modswitch = compute_modswitch(s2_total, pow(2, modswitch_pbits), pow(2, 50), h)
    print("log p = ", math.log2(compute_exp(modswitch_pbits, math.log2(p), s2_modswitch)))

# 8 GB
def test_probability_bsgs_8GB():
    N = 4096
    p = 65537 # 40961 # 65537
    h = 2731
    # h = 256
    # h = 128

    ell = 3
    Bg = 0x01 << 20
    PP = 0x01 << 20

    N1 = 128 # 64

    ell_gs = 2 # 3
    Bg_gs = 0x01 << 18 # 0x01 << 14
    PP_gs = 0x01 << 20 # 0x01 << 14
    # ell_gs = 49 # 3
    # Bg_gs = 0x01 << 1 # 0x01 << 12
    # PP_gs = 0x01 << 1 # 0x01 << 14

    sigma2_first_dim = compute_sigma2_bsgs(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)
    # sigma2_first_dim = compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(52, math.log2(p), sigma2_first_dim)))

    # if we add a modswitch here
    modswitch_p1bits = 50
    sigma2_first_dim = compute_modswitch(sigma2_first_dim, pow(2, modswitch_p1bits), pow(2, 56), h)
    # print("log p = ", math.log2(compute_exp(modswitch_p1bits, math.log2(p), sigma2_first_dim)))

    # external product
    l_exP = 2
    Bg_exP = pow(2, 16) # pow(2, 18)
    base_exP = pow(2, 18) # pow(2, 20)
    sigma2_second_dim = compute_s2_second_dimension(N, sigma2_first_dim, l_exP, Bg_exP, base_exP)


    # packing lwes
    l_atm = 4 # 3 # 4
    Bg_atm = pow(2, 13) # pow(2, 14) # pow(2, 10)
    base_atm = 0 # pow(2, 14) # pow(2, 16)
    s2_total = compute_packing_lwes_formula_new(16 * 4 * 2 * 4, N, l_atm, Bg_atm, base_atm, sigma2_second_dim)

    print("s2_total = ", s2_total)
    print("log p = ", math.log2(compute_exp(47, math.log2(p), s2_total)))

    modswitch_pbits = 26
    s2_modswitch = compute_modswitch(s2_total, pow(2, modswitch_pbits), pow(2, 50), h)
    print("log p = ", math.log2(compute_exp(modswitch_pbits, math.log2(p), s2_modswitch)))

def test_probability_bsgs_stateless_1GB():
    N = 4096
    p = 65537 # 40961 # 65537
    h = 2731
    # h = 256
    # h = 128

    ell = 3
    Bg = 0x01 << 20
    PP = 0x01 << 20

    N1 = 512 # 128 # 64

    ell_gs = 2 # 3
    Bg_gs = 0x01 << 18 # 0x01 << 14
    PP_gs = 0x01 << 20 # 0x01 << 14
    # ell_gs = 49 # 3
    # Bg_gs = 0x01 << 1 # 0x01 << 12
    # PP_gs = 0x01 << 1 # 0x01 << 14

    sigma2_first_dim = compute_sigma2_bsgs_stateless(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)
    # sigma2_first_dim = compute_sigma2_bsgs_offlineKey(N, p, ell, Bg, PP, ell_gs, Bg_gs, PP_gs, N1, h)

    print("sigma2_first_dim = ", sigma2_first_dim)
    print("log p = ", math.log2(compute_exp(52, math.log2(p), sigma2_first_dim)))

    # if we add a modswitch here
    modswitch_p1bits = 50
    sigma2_first_dim = compute_modswitch(sigma2_first_dim, pow(2, modswitch_p1bits), pow(2, 56), h)
    # print("log p = ", math.log2(compute_exp(modswitch_p1bits, math.log2(p), sigma2_first_dim)))

    # external product
    l_exP = 2
    Bg_exP = pow(2, 16) # pow(2, 18)
    base_exP = pow(2, 18) # pow(2, 20)
    sigma2_second_dim = compute_s2_second_dimension(N, sigma2_first_dim, l_exP, Bg_exP, base_exP)


    # packing lwes
    l_atm = 4 # 3 # 4
    Bg_atm = pow(2, 13) # pow(2, 14) # pow(2, 10)
    base_atm = 0 # pow(2, 14) # pow(2, 16)
    l_atm = 2 # 3
    Bg_atm = pow(2, 16) # pow(2, 14)
    base_atm = pow(2, 18) # pow(2, 14)
    s2_total = compute_packing_lwes_formula_new(16 * 4, N, l_atm, Bg_atm, base_atm, sigma2_second_dim)

    print("s2_total = ", s2_total)
    print("log p = ", math.log2(compute_exp(47, math.log2(p), s2_total)))

    modswitch_pbits = 26
    s2_modswitch = compute_modswitch(s2_total, pow(2, modswitch_pbits), pow(2, 50), h)
    print("log p = ", math.log2(compute_exp(modswitch_pbits, math.log2(p), s2_modswitch)))


if __name__ == '__main__':

    # test_probability()
    # test_probability_bsgs_36MB()
    # test_probability_bsgs() # 256 MB
    test_probability_bsgs_1GB()
    # test_probability_bsgs_2GB()
    # test_probability_bsgs_8GB()