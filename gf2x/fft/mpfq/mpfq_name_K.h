#ifndef MPFQ_NAME_K_H_
#define MPFQ_NAME_K_H_

/* Automatically generated file.
 *
 * This header file is just wrap-up code for accessing a global finite
 * field with convenient aliases.
 *
 * Note that this file is automatically generated from the mpfq api, and
 * is therefore guaranteed to contain all the api functions usable in the
 * current api_extensions.
 */

#ifndef MPFQ_LAST_GENERATED_TAG
#error "Please include an mpfq-generated header first"
#endif

/* cpp has its infelicities. Yes the extra step is needed */
#ifndef MPFQ_CONCAT4
#define MPFQ_CONCAT4(X,Y,Z,T) X ## Y ## Z ## T
#endif

#ifndef MPFQ_CREATE_FUNCTION_NAME
#define MPFQ_CREATE_FUNCTION_NAME(TAG,NAME) MPFQ_CONCAT4(mpfq_,TAG,_,NAME)
#endif

#define Kcommon_name_(X) MPFQ_CREATE_FUNCTION_NAME(MPFQ_LAST_GENERATED_TAG,X)

#define Kfield	Kcommon_name_(field)
#define Kdst_field	Kcommon_name_(dst_field)
#define Ksrc_field	Kcommon_name_(src_field)

#define Kelt	Kcommon_name_(elt)
#define Kdst_elt	Kcommon_name_(dst_elt)
#define Ksrc_elt	Kcommon_name_(src_elt)

#define Kelt_ur	Kcommon_name_(elt_ur)
#define Kdst_elt_ur	Kcommon_name_(dst_elt_ur)
#define Ksrc_elt_ur	Kcommon_name_(src_elt_ur)

#define Kvec	Kcommon_name_(vec)
#define Kdst_vec	Kcommon_name_(dst_vec)
#define Ksrc_vec	Kcommon_name_(src_vec)

#define Kvec_ur	Kcommon_name_(vec_ur)
#define Kdst_vec_ur	Kcommon_name_(dst_vec_ur)
#define Ksrc_vec_ur	Kcommon_name_(src_vec_ur)

#define Kpoly	Kcommon_name_(poly)
#define Kdst_poly	Kcommon_name_(dst_poly)
#define Ksrc_poly	Kcommon_name_(src_poly)






#define Kimpl_name()	Kcommon_name_(impl_name) ()
#define Kimpl_max_characteristic_bits()	Kcommon_name_(impl_max_characteristic_bits) ()
#define Kimpl_max_degree()	Kcommon_name_(impl_max_degree) ()



#define Kfield_characteristic(a0)	Kcommon_name_(field_characteristic) (K,a0)
#define Kfield_characteristic_srcptr()	Kcommon_name_(field_characteristic_srcptr) (K)
#define Kfield_characteristic_bits()	Kcommon_name_(field_characteristic_bits) (K)
#define Kfield_degree()	Kcommon_name_(field_degree) (K)

#define Kfield_init()	Kcommon_name_(field_init) (K)
#define Kfield_clear()	Kcommon_name_(field_clear) (K)
#define Kfield_specify(a0,a1)	Kcommon_name_(field_specify) (K,a0,a1)
#define Kfield_setopt(a0,a1)	Kcommon_name_(field_setopt) (K,a0,a1)

#define Kinit(a0)	Kcommon_name_(init) (K,a0)
#define Kclear(a0)	Kcommon_name_(clear) (K,a0)
#define Kelt_stride()	Kcommon_name_(elt_stride) (K)

#define Kset(a0,a1)	Kcommon_name_(set) (K,a0,a1)
#define Kset_ui(a0,a1)	Kcommon_name_(set_ui) (K,a0,a1)
#define Kset_zero(a0)	Kcommon_name_(set_zero) (K,a0)
#define Kget_ui(a0)	Kcommon_name_(get_ui) (K,a0)
#define Kset_mpn(a0,a1,a2)	Kcommon_name_(set_mpn) (K,a0,a1,a2)
#define Kset_mpz(a0,a1)	Kcommon_name_(set_mpz) (K,a0,a1)
#define Kget_mpn(a0,a1)	Kcommon_name_(get_mpn) (K,a0,a1)
#define Kget_mpz(a0,a1)	Kcommon_name_(get_mpz) (K,a0,a1)

#define Kset_uipoly(a0,a1)	Kcommon_name_(set_uipoly) (K,a0,a1)
#define Kset_uipoly_wide(a0,a1,a2)	Kcommon_name_(set_uipoly_wide) (K,a0,a1,a2)
#define Kget_uipoly(a0)	Kcommon_name_(get_uipoly) (K,a0)
#define Kget_uipoly_wide(a0,a1)	Kcommon_name_(get_uipoly_wide) (K,a0,a1)


#define Krandom(a0,a1)	Kcommon_name_(random) (K,a0,a1)
#define Krandom2(a0,a1)	Kcommon_name_(random2) (K,a0,a1)

#define Kadd(a0,a1,a2)	Kcommon_name_(add) (K,a0,a1,a2)
#define Ksub(a0,a1,a2)	Kcommon_name_(sub) (K,a0,a1,a2)
#define Kneg(a0,a1)	Kcommon_name_(neg) (K,a0,a1)
#define Kmul(a0,a1,a2)	Kcommon_name_(mul) (K,a0,a1,a2)
#define Ksqr(a0,a1)	Kcommon_name_(sqr) (K,a0,a1)
#define Kis_sqr(a0)	Kcommon_name_(is_sqr) (K,a0)
#define Ksqrt(a0,a1)	Kcommon_name_(sqrt) (K,a0,a1)
#define Kpow(a0,a1,a2,a3)	Kcommon_name_(pow) (K,a0,a1,a2,a3)
#define Kpowz(a0,a1,a2)	Kcommon_name_(powz) (K,a0,a1,a2)
#define Kfrobenius(a0,a1)	Kcommon_name_(frobenius) (K,a0,a1)
#define Kadd_ui(a0,a1,a2)	Kcommon_name_(add_ui) (K,a0,a1,a2)
#define Ksub_ui(a0,a1,a2)	Kcommon_name_(sub_ui) (K,a0,a1,a2)
#define Kmul_ui(a0,a1,a2)	Kcommon_name_(mul_ui) (K,a0,a1,a2)
#define Knormalize(a0)	Kcommon_name_(normalize) (K,a0)
#define Kadd_uipoly(a0,a1,a2)	Kcommon_name_(add_uipoly) (K,a0,a1,a2)
#define Ksub_uipoly(a0,a1,a2)	Kcommon_name_(sub_uipoly) (K,a0,a1,a2)
#define Kmul_uipoly(a0,a1,a2)	Kcommon_name_(mul_uipoly) (K,a0,a1,a2)
#define Kinv(a0,a1)	Kcommon_name_(inv) (K,a0,a1)
#define Kas_solve(a0,a1)	Kcommon_name_(as_solve) (K,a0,a1)
#define Ktrace(a0)	Kcommon_name_(trace) (K,a0)
#define Khadamard(a0,a1,a2,a3)	Kcommon_name_(hadamard) (K,a0,a1,a2,a3)


#define Kelt_ur_init(a0)	Kcommon_name_(elt_ur_init) (K,a0)
#define Kelt_ur_clear(a0)	Kcommon_name_(elt_ur_clear) (K,a0)
#define Kelt_ur_stride()	Kcommon_name_(elt_ur_stride) (K)
#define Kelt_ur_set(a0,a1)	Kcommon_name_(elt_ur_set) (K,a0,a1)
#define Kelt_ur_set_elt(a0,a1)	Kcommon_name_(elt_ur_set_elt) (K,a0,a1)
#define Kelt_ur_set_zero(a0)	Kcommon_name_(elt_ur_set_zero) (K,a0)
#define Kelt_ur_set_ui(a0,a1)	Kcommon_name_(elt_ur_set_ui) (K,a0,a1)
#define Kelt_ur_add(a0,a1,a2)	Kcommon_name_(elt_ur_add) (K,a0,a1,a2)
#define Kelt_ur_neg(a0,a1)	Kcommon_name_(elt_ur_neg) (K,a0,a1)
#define Kelt_ur_sub(a0,a1,a2)	Kcommon_name_(elt_ur_sub) (K,a0,a1,a2)
#define Kmul_ur(a0,a1,a2)	Kcommon_name_(mul_ur) (K,a0,a1,a2)
#define Ksqr_ur(a0,a1)	Kcommon_name_(sqr_ur) (K,a0,a1)
#define Kreduce(a0,a1)	Kcommon_name_(reduce) (K,a0,a1)
#define Kaddmul_si_ur(a0,a1,a2)	Kcommon_name_(addmul_si_ur) (K,a0,a1,a2)


#define Kcmp(a0,a1)	Kcommon_name_(cmp) (K,a0,a1)
#define Kcmp_ui(a0,a1)	Kcommon_name_(cmp_ui) (K,a0,a1)
#define Kis_zero(a0)	Kcommon_name_(is_zero) (K,a0)


#define Kmgy_enc(a0,a1)	Kcommon_name_(mgy_enc) (K,a0,a1)
#define Kmgy_dec(a0,a1)	Kcommon_name_(mgy_dec) (K,a0,a1)


#define Kasprint(a0,a1)	Kcommon_name_(asprint) (K,a0,a1)
#define Kfprint(a0,a1)	Kcommon_name_(fprint) (K,a0,a1)
#define Kcxx_out(a0,a1)	Kcommon_name_(cxx_out) (K,a0,a1)
#define Kprint(a0)	Kcommon_name_(print) (K,a0)
#define Ksscan(a0,a1)	Kcommon_name_(sscan) (K,a0,a1)
#define Kfscan(a0,a1)	Kcommon_name_(fscan) (K,a0,a1)
#define Kcxx_in(a0,a1)	Kcommon_name_(cxx_in) (K,a0,a1)
#define Kscan(a0)	Kcommon_name_(scan) (K,a0)
#define Kread(a0,a1)	Kcommon_name_(read) (K,a0,a1)
#define Kimportdata(a0,a1,a2,a3)	Kcommon_name_(importdata) (K,a0,a1,a2,a3)
#define Kwrite(a0,a1)	Kcommon_name_(write) (K,a0,a1)
#define Kexportdata(a0,a1,a2,a3)	Kcommon_name_(exportdata) (K,a0,a1,a2,a3)





#define Kvec_init(a0,a1)	Kcommon_name_(vec_init) (K,a0,a1)
#define Kvec_reinit(a0,a1,a2)	Kcommon_name_(vec_reinit) (K,a0,a1,a2)
#define Kvec_clear(a0,a1)	Kcommon_name_(vec_clear) (K,a0,a1)
#define Kvec_set(a0,a1,a2)	Kcommon_name_(vec_set) (K,a0,a1,a2)
#define Kvec_set_zero(a0,a1)	Kcommon_name_(vec_set_zero) (K,a0,a1)
#define Kvec_setcoeff(a0,a1,a2)	Kcommon_name_(vec_setcoeff) (K,a0,a1,a2)
#define Kvec_setcoeff_ui(a0,a1,a2)	Kcommon_name_(vec_setcoeff_ui) (K,a0,a1,a2)
#define Kvec_getcoeff(a0,a1,a2)	Kcommon_name_(vec_getcoeff) (K,a0,a1,a2)
#define Kvec_add(a0,a1,a2,a3)	Kcommon_name_(vec_add) (K,a0,a1,a2,a3)
#define Kvec_neg(a0,a1,a2)	Kcommon_name_(vec_neg) (K,a0,a1,a2)
#define Kvec_rev(a0,a1,a2)	Kcommon_name_(vec_rev) (K,a0,a1,a2)
#define Kvec_sub(a0,a1,a2,a3)	Kcommon_name_(vec_sub) (K,a0,a1,a2,a3)
#define Kvec_scal_mul(a0,a1,a2,a3)	Kcommon_name_(vec_scal_mul) (K,a0,a1,a2,a3)
#define Kvec_conv(a0,a1,a2,a3,a4)	Kcommon_name_(vec_conv) (K,a0,a1,a2,a3,a4)
#define Kvec_random(a0,a1,a2)	Kcommon_name_(vec_random) (K,a0,a1,a2)
#define Kvec_random2(a0,a1,a2)	Kcommon_name_(vec_random2) (K,a0,a1,a2)
#define Kvec_cmp(a0,a1,a2)	Kcommon_name_(vec_cmp) (K,a0,a1,a2)
#define Kvec_is_zero(a0,a1)	Kcommon_name_(vec_is_zero) (K,a0,a1)
#define Kvec_subvec(a0,a1)	Kcommon_name_(vec_subvec) (K,a0,a1)
#define Kvec_subvec_const(a0,a1)	Kcommon_name_(vec_subvec_const) (K,a0,a1)
#define Kvec_coeff_ptr(a0,a1)	Kcommon_name_(vec_coeff_ptr) (K,a0,a1)
#define Kvec_coeff_ptr_const(a0,a1)	Kcommon_name_(vec_coeff_ptr_const) (K,a0,a1)
#define Kvec_asprint(a0,a1,a2)	Kcommon_name_(vec_asprint) (K,a0,a1,a2)
#define Kvec_fprint(a0,a1,a2)	Kcommon_name_(vec_fprint) (K,a0,a1,a2)
#define Kvec_print(a0,a1)	Kcommon_name_(vec_print) (K,a0,a1)
#define Kvec_sscan(a0,a1,a2)	Kcommon_name_(vec_sscan) (K,a0,a1,a2)
#define Kvec_fscan(a0,a1,a2)	Kcommon_name_(vec_fscan) (K,a0,a1,a2)
#define Kvec_scan(a0,a1)	Kcommon_name_(vec_scan) (K,a0,a1)
#define Kvec_cxx_out(a0,a1,a2)	Kcommon_name_(vec_cxx_out) (K,a0,a1,a2)
#define Kvec_cxx_in(a0,a1,a2)	Kcommon_name_(vec_cxx_in) (K,a0,a1,a2)
#define Kvec_read(a0,a1,a2)	Kcommon_name_(vec_read) (K,a0,a1,a2)
#define Kvec_write(a0,a1,a2)	Kcommon_name_(vec_write) (K,a0,a1,a2)
#define Kvec_import(a0,a1,a2,a3)	Kcommon_name_(vec_import) (K,a0,a1,a2,a3)
#define Kvec_export(a0,a1,a2,a3)	Kcommon_name_(vec_export) (K,a0,a1,a2,a3)
#define Kvec_hamming_weight(a0,a1)	Kcommon_name_(vec_hamming_weight) (K,a0,a1)
#define Kvec_find_first_set(a0,a1)	Kcommon_name_(vec_find_first_set) (K,a0,a1)

#define Kvec_simd_hamming_weight(a0,a1)	Kcommon_name_(vec_simd_hamming_weight) (K,a0,a1)
#define Kvec_simd_find_first_set(a0,a1)	Kcommon_name_(vec_simd_find_first_set) (K,a0,a1)


#define Kvec_ur_init(a0,a1)	Kcommon_name_(vec_ur_init) (K,a0,a1)
#define Kvec_ur_set_zero(a0,a1)	Kcommon_name_(vec_ur_set_zero) (K,a0,a1)
#define Kvec_ur_set_vec(a0,a1,a2)	Kcommon_name_(vec_ur_set_vec) (K,a0,a1,a2)
#define Kvec_ur_reinit(a0,a1,a2)	Kcommon_name_(vec_ur_reinit) (K,a0,a1,a2)
#define Kvec_ur_clear(a0,a1)	Kcommon_name_(vec_ur_clear) (K,a0,a1)
#define Kvec_ur_set(a0,a1,a2)	Kcommon_name_(vec_ur_set) (K,a0,a1,a2)
#define Kvec_ur_setcoeff(a0,a1,a2)	Kcommon_name_(vec_ur_setcoeff) (K,a0,a1,a2)
#define Kvec_ur_getcoeff(a0,a1,a2)	Kcommon_name_(vec_ur_getcoeff) (K,a0,a1,a2)
#define Kvec_ur_add(a0,a1,a2,a3)	Kcommon_name_(vec_ur_add) (K,a0,a1,a2,a3)
#define Kvec_ur_sub(a0,a1,a2,a3)	Kcommon_name_(vec_ur_sub) (K,a0,a1,a2,a3)
#define Kvec_ur_neg(a0,a1,a2)	Kcommon_name_(vec_ur_neg) (K,a0,a1,a2)
#define Kvec_ur_rev(a0,a1,a2)	Kcommon_name_(vec_ur_rev) (K,a0,a1,a2)
#define Kvec_scal_mul_ur(a0,a1,a2,a3)	Kcommon_name_(vec_scal_mul_ur) (K,a0,a1,a2,a3)
#define Kvec_conv_ur(a0,a1,a2,a3,a4)	Kcommon_name_(vec_conv_ur) (K,a0,a1,a2,a3,a4)
#define Kvec_reduce(a0,a1,a2)	Kcommon_name_(vec_reduce) (K,a0,a1,a2)
#define Kvec_ur_subvec(a0,a1)	Kcommon_name_(vec_ur_subvec) (K,a0,a1)
#define Kvec_ur_subvec_const(a0,a1)	Kcommon_name_(vec_ur_subvec_const) (K,a0,a1)
#define Kvec_ur_coeff_ptr(a0,a1)	Kcommon_name_(vec_ur_coeff_ptr) (K,a0,a1)
#define Kvec_ur_coeff_ptr_const(a0,a1)	Kcommon_name_(vec_ur_coeff_ptr_const) (K,a0,a1)

#define Kvec_elt_stride(a0)	Kcommon_name_(vec_elt_stride) (K,a0)

#define Kvec_ur_elt_stride(a0)	Kcommon_name_(vec_ur_elt_stride) (K,a0)




#define Kpoly_init(a0,a1)	Kcommon_name_(poly_init) (K,a0,a1)
#define Kpoly_clear(a0)	Kcommon_name_(poly_clear) (K,a0)
#define Kpoly_set(a0,a1)	Kcommon_name_(poly_set) (K,a0,a1)
#define Kpoly_setmonic(a0,a1)	Kcommon_name_(poly_setmonic) (K,a0,a1)
#define Kpoly_setcoeff(a0,a1,a2)	Kcommon_name_(poly_setcoeff) (K,a0,a1,a2)
#define Kpoly_setcoeff_ui(a0,a1,a2)	Kcommon_name_(poly_setcoeff_ui) (K,a0,a1,a2)
#define Kpoly_getcoeff(a0,a1,a2)	Kcommon_name_(poly_getcoeff) (K,a0,a1,a2)
#define Kpoly_deg(a0)	Kcommon_name_(poly_deg) (K,a0)
#define Kpoly_add(a0,a1,a2)	Kcommon_name_(poly_add) (K,a0,a1,a2)
#define Kpoly_sub(a0,a1,a2)	Kcommon_name_(poly_sub) (K,a0,a1,a2)
#define Kpoly_set_ui(a0,a1)	Kcommon_name_(poly_set_ui) (K,a0,a1)
#define Kpoly_add_ui(a0,a1,a2)	Kcommon_name_(poly_add_ui) (K,a0,a1,a2)
#define Kpoly_sub_ui(a0,a1,a2)	Kcommon_name_(poly_sub_ui) (K,a0,a1,a2)
#define Kpoly_neg(a0,a1)	Kcommon_name_(poly_neg) (K,a0,a1)
#define Kpoly_scal_mul(a0,a1,a2)	Kcommon_name_(poly_scal_mul) (K,a0,a1,a2)
#define Kpoly_mul(a0,a1,a2)	Kcommon_name_(poly_mul) (K,a0,a1,a2)
#define Kpoly_divmod(a0,a1,a2,a3)	Kcommon_name_(poly_divmod) (K,a0,a1,a2,a3)
#define Kpoly_precomp_mod(a0,a1)	Kcommon_name_(poly_precomp_mod) (K,a0,a1)
#define Kpoly_mod_pre(a0,a1,a2,a3)	Kcommon_name_(poly_mod_pre) (K,a0,a1,a2,a3)
#define Kpoly_gcd(a0,a1,a2)	Kcommon_name_(poly_gcd) (K,a0,a1,a2)
#define Kpoly_xgcd(a0,a1,a2,a3,a4)	Kcommon_name_(poly_xgcd) (K,a0,a1,a2,a3,a4)
#define Kpoly_random(a0,a1,a2)	Kcommon_name_(poly_random) (K,a0,a1,a2)
#define Kpoly_random2(a0,a1,a2)	Kcommon_name_(poly_random2) (K,a0,a1,a2)
#define Kpoly_cmp(a0,a1)	Kcommon_name_(poly_cmp) (K,a0,a1)
#define Kpoly_asprint(a0,a1)	Kcommon_name_(poly_asprint) (K,a0,a1)
#define Kpoly_fprint(a0,a1)	Kcommon_name_(poly_fprint) (K,a0,a1)
#define Kpoly_print(a0)	Kcommon_name_(poly_print) (K,a0)
#define Kpoly_sscan(a0,a1)	Kcommon_name_(poly_sscan) (K,a0,a1)
#define Kpoly_fscan(a0,a1)	Kcommon_name_(poly_fscan) (K,a0,a1)
#define Kpoly_scan(a0)	Kcommon_name_(poly_scan) (K,a0)
#define Kpoly_cxx_out(a0,a1)	Kcommon_name_(poly_cxx_out) (K,a0,a1)
#define Kpoly_cxx_in(a0,a1)	Kcommon_name_(poly_cxx_in) (K,a0,a1)



#define Ksimd_groupsize()	Kcommon_name_(simd_groupsize) (K)
#define Ksimd_hamming_weight(a0)	Kcommon_name_(simd_hamming_weight) (K,a0)
#define Ksimd_find_first_set(a0)	Kcommon_name_(simd_find_first_set) (K,a0)
#define Ksimd_get_ui_at(a0,a1)	Kcommon_name_(simd_get_ui_at) (K,a0,a1)
#define Ksimd_set_ui_at(a0,a1,a2)	Kcommon_name_(simd_set_ui_at) (K,a0,a1,a2)
#define Ksimd_add_ui_at(a0,a1,a2,a3)	Kcommon_name_(simd_add_ui_at) (K,a0,a1,a2,a3)
#define Ksimd_set_ui_all(a0,a1)	Kcommon_name_(simd_set_ui_all) (K,a0,a1)
#define Kadd_dotprod(a0,a1,a2,a3)	Kcommon_name_(add_dotprod) (K,a0,a1,a2,a3)
#define Kmul_constant_ui(a0,a1,a2)	Kcommon_name_(mul_constant_ui) (K,a0,a1,a2)


#define Kmember_template_add_dotprod(a0,a1,a2,a3,a4)	Kcommon_name_(member_template_add_dotprod) (K,a0,a1,a2,a3,a4)
#define Kmember_template_addmul_tiny(a0,a1,a2,a3,a4)	Kcommon_name_(member_template_addmul_tiny) (K,a0,a1,a2,a3,a4)
#define Kmember_template_transpose(a0,a1,a2)	Kcommon_name_(member_template_transpose) (K,a0,a1,a2)





#define Koo_field_init()	Kcommon_name_(oo_field_init) (K)
#define Koo_field_clear()	Kcommon_name_(oo_field_clear) (K)


/* customary link reference to the field -- forces good habit of defining
   it somewhere */
/* another customary shorthand */
#define	Kdegree	Kfield_degree()


#endif  /* MPFQ_NAME_K_H_ */
