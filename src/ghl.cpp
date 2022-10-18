#include "test.h"
#include "bench.h"
#include "common.h"
#include <sys/random.h>

void ghl_sample_message(params::poly_2q & m) {
	// Sample a short polynomial.
	std::array < mpz_t, 2 * DEGREE > coeffs;
	uint64_t buf;

	size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], bits_in_moduli_product << 2);
	}

	for (size_t j = 0; j < params::poly_2q::degree / 2; j += 32) {
		getrandom(&buf, sizeof(buf), 0);
		for (size_t k = 0; k < 64; k += 2) {
			mpz_set_ui(coeffs[j + k / 2], (buf >> k) % PRIMEP);
		}
	}
	m.mpz2poly(coeffs);

	for (size_t i = 0; i < params::poly_p::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

void ghl_encode(params::poly_2q & m) {
	// Sample a short polynomial.
	mpz_t q;
	std::array < mpz_t, 2 * DEGREE > coeffs, encoding;

	mpz_init(q);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_2q::bits_in_moduli_product() << 2);
		mpz_init2(encoding[i], params::poly_2q::bits_in_moduli_product() << 2);
	}

	m.poly2mpz(coeffs);
	/* Encode for GHL encryption. */
	mpz_set_str(q, PRIMEQ, 10);
	for (size_t j = 0; j < params::poly_2q::degree; j += 2) {
		mpz_set(encoding[j + 1], coeffs[j / 2]);
		mpz_mul_ui(encoding[j], coeffs[j / 2], DELTA);
		mpz_mod(encoding[j], encoding[j], q);
	}
	m.mpz2poly(encoding);

	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
	mpz_clear(q);
}

void ghl_decode(params::poly_2q & m) {
	// Sample a short polynomial.
	mpz_t q, qDivBy2;
	std::array < mpz_t, 2 * DEGREE > coeffs, encoding;

	mpz_init(q);
	mpz_init(qDivBy2);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_2q::bits_in_moduli_product() << 2);
		mpz_init2(encoding[i], params::poly_2q::bits_in_moduli_product() << 2);
	}
	mpz_fdiv_q_2exp(qDivBy2, params::poly_2q::moduli_product(), 1);

	m.poly2mpz(encoding);
	/* Decode for GHL encryption. */
	mpz_set_str(q, PRIMEQ, 10);
	for (size_t j = 0; j < params::poly_2q::degree; j += 2) {
		mpz_set(coeffs[j / 2], encoding[j]);
		mpz_mul_ui(encoding[j], encoding[j], DELTA);
		mpz_sub(encoding[j], encoding[j + 1], encoding[j]);
		mpz_mod(encoding[j], encoding[j], q);
		mpz_mod_ui(encoding[j], encoding[j], DELTA);
		mpz_sub(coeffs[j / 2], coeffs[j / 2], encoding[j]);
		// TODO: make encoding more general by solving rounding below
		mpz_div_ui(coeffs[j / 2], coeffs[j / 2], DELTA_R);
		//gmp_printf("coeff %d: %Zd \n", j, coeffs[j/2]);
	}
	m.mpz2poly(coeffs);

	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
	mpz_clear(q);
	mpz_clear(qDivBy2);
}

// Generate a key pair.
void ghl_keygen(ghlkey_t & pk, params::poly_2q & sk) {
	params::poly_2q e = nfl::ZO_dist();
	pk.a = nfl::uniform();

	sk = nfl::ZO_dist();
	sk.ntt_pow_phi();
	e.ntt_pow_phi();

	pk.b = pk.a * sk + e;
}

void ghl_keyshare(params::poly_2q s[], size_t shares, params::poly_2q & sk) {
	params::poly_2q t = sk;
	for (size_t i = 1; i < shares; i++) {
		s[i] = nfl::uniform();
		s[i].ntt_pow_phi();
		t = t - s[i];
	}
	s[0] = t;
}

void ghl_encrypt(ghlenc_t & c, ghlkey_t & pk, params::poly_2q & m) {
	std::array < mpz_t, 2 * DEGREE > coeffs;
	params::poly_2q e1 = nfl::ZO_dist();
	params::poly_2q e2 = nfl::ZO_dist();
	params::poly_2q r = nfl::ZO_dist();

	e1.ntt_pow_phi();
	e2.ntt_pow_phi();
	r.ntt_pow_phi();

	c.u = pk.a * r + e1;
	c.v = pk.b * r + e2;

	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_2q::bits_in_moduli_product() << 2);
	}

	m.poly2mpz(coeffs);
	r.mpz2poly(coeffs);
	r.ntt_pow_phi();
	c.v = c.v + r;

	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

void ghl_decrypt(params::poly_2q & m, ghlenc_t & c, params::poly_2q & sk) {
	std::array < mpz_t, 2 * DEGREE > coeffs;
	params::poly_2q t = c.v - sk * c.u;
	mpz_t qDivBy2;

	mpz_init(qDivBy2);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_2q::bits_in_moduli_product() << 2);
	}
	mpz_fdiv_q_2exp(qDivBy2, params::poly_2q::moduli_product(), 1);

	// Reduce the coefficients
	t.invntt_pow_invphi();
	t.poly2mpz(coeffs);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		util::center(coeffs[i], coeffs[i], params::poly_2q::moduli_product(),
				qDivBy2);
	}
	m.mpz2poly(coeffs);

	mpz_clear(qDivBy2);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

void ghl_distdec(params::poly_2q & tj, ghlenc_t & c, params::poly_2q & sj) {
	std::array < mpz_t, 2 * DEGREE > coeffs;
	params::poly_2q mj, Ej;
	mpz_t qDivBy2, bound;

	mpz_init(qDivBy2);
	mpz_init(bound);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_2q::bits_in_moduli_product() << 2);
	}
	mpz_fdiv_q_2exp(qDivBy2, params::poly_2q::moduli_product(), 1);
	mpz_set_str(bound, BOUND_D, 10);

	Ej = nfl::uniform();
	Ej.poly2mpz(coeffs);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		util::center(coeffs[i], coeffs[i], params::poly_2q::moduli_product(),
				qDivBy2);
		mpz_mod(coeffs[i], coeffs[i], bound);
	}
	Ej.mpz2poly(coeffs);
	Ej.ntt_pow_phi();
	mj = sj * c.u;
	tj = mj + Ej;
}

void ghl_comb(params::poly_2q & m, ghlenc_t & c, params::poly_2q t[],
		size_t shares) {
	std::array < mpz_t, 2 * DEGREE > coeffs;
	params::poly_2q v;
	mpz_t qDivBy2;

	mpz_init(qDivBy2);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_init2(coeffs[i], params::poly_2q::bits_in_moduli_product() << 2);
	}
	mpz_fdiv_q_2exp(qDivBy2, params::poly_2q::moduli_product(), 1);

	v = c.v - t[0];
	for (size_t i = 1; i < shares; i++) {
		v = v - t[i];
	}
	v.invntt_pow_invphi();
	// Reduce the coefficients modulo p
	v.poly2mpz(coeffs);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		util::center(coeffs[i], coeffs[i], params::poly_2q::moduli_product(),
				qDivBy2);
		mpz_mod_ui(coeffs[i], coeffs[i], PRIMEP);
	}
	m.mpz2poly(coeffs);

	mpz_clear(qDivBy2);
	for (size_t i = 0; i < params::poly_2q::degree; i++) {
		mpz_clear(coeffs[i]);
	}
}

#ifdef MAIN
static void test() {
	ghlkey_t pk;
	params::poly_2q sk, s[PARTIES], t[PARTIES], acc, m, _m;
	ghlenc_t c;

	ghl_keygen(pk, sk);
	ghl_sample_message(m);

	TEST_BEGIN("GHL encoding is consistent") {
		_m = m;
		ghl_encode(_m);
		ghl_decode(_m);
		TEST_ASSERT(m - _m == 0, end);
	} TEST_END;

  end:
	return;
}

static void bench() {
	ghlkey_t pk;
	params::poly_2q sk, s[PARTIES], t[PARTIES], acc, m, _m;
	ghlenc_t c;

	ghl_keygen(pk, sk);

	BENCH_BEGIN("ghl_sample_message") {
		BENCH_ADD(ghl_sample_message(m));
	} BENCH_END;

	BENCH_BEGIN("ghl_encrypt") {
		BENCH_ADD(ghl_encrypt(c, pk, m));
	} BENCH_END;

	BENCH_BEGIN("ghl_decrypt") {
		ghl_encrypt(c, pk, m);
		BENCH_ADD(ghl_decrypt(_m, c, sk));
	} BENCH_END;

	ghl_keyshare(t, PARTIES, sk);
	BENCH_BEGIN("ghl_distdec") {
		ghl_encrypt(c, pk, m);
		BENCH_ADD(ghl_distdec(t[0], c, s[0]));
	} BENCH_END;

	for (size_t i = 1; i < PARTIES; i++) {
		ghl_distdec(t[i], c, s[i]);
	}

	BENCH_BEGIN("ghl_comb") {
		BENCH_ADD(ghl_comb(_m, c, t, PARTIES));
	} BENCH_END;
}

int main(int argc, char *arv[]) {
	printf("\n** Tests for GHL encryption:\n\n");
	test();

	printf("\n** Benchmarks for GHL encryption:\n\n");
	bench();

	return 0;
}
#endif
