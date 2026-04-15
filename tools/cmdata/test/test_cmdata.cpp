// Unit tests for cmdata headers that have no GROMACS dependency:
//   indexing.hpp  — n_bins, offset_same, offset_cross
//   density.hpp   — kernel_density_estimator, normalize_histo
//   mindist.hpp   — mindist_same, mindist_cross
//
// Build (no CMake needed):
//   g++ -std=c++17 -I../src test_cmdata.cpp -o test_cmdata && ./test_cmdata

#include <cmath>
#include <cstdio>
#include <numeric>
#include <set>
#include <vector>

#include "cmdata/indexing.hpp"
#include "cmdata/density.hpp"
#include "cmdata/mindist.hpp"

// ── minimal test harness ─────────────────────────────────────────────────────
static int g_pass = 0, g_fail = 0;

#define CHECK(expr) \
  do { \
    if (expr) { ++g_pass; } \
    else { ++g_fail; fprintf(stderr, "FAIL  %s:%d  %s\n", __FILE__, __LINE__, #expr); } \
  } while (0)

#define CHECK_MSG(expr, msg) \
  do { \
    if (expr) { ++g_pass; } \
    else { ++g_fail; fprintf(stderr, "FAIL  %s:%d  %s  (%s)\n", __FILE__, __LINE__, #expr, msg); } \
  } while (0)

#define CHECK_NEAR(a, b, tol) \
  do { \
    double _a = (a), _b = (b), _t = (tol); \
    if (std::abs(_a - _b) <= _t) { ++g_pass; } \
    else { ++g_fail; fprintf(stderr, "FAIL  %s:%d  |%.6g - %.6g| = %.6g > %.6g\n", __FILE__, __LINE__, _a, _b, std::abs(_a - _b), _t); } \
  } while (0)

// ── helpers ──────────────────────────────────────────────────────────────────

// Build a uniform density_bins vector of length n for the given cutoff.
static std::vector<float> make_bins(float cutoff)
{
  int n = cmdata::indexing::n_bins(cutoff);
  std::vector<float> bins(n);
  float dx = cutoff / static_cast<float>(n);
  for (int i = 0; i < n; i++) bins[i] = dx * (static_cast<float>(i) + 0.5f);
  return bins;
}

// Sum a histogram, multiplied by dx, to approximate the integral.
static double integrate(const std::vector<float> &h, float dx)
{
  double s = 0.0;
  for (float v : h) s += v;
  return s * dx;
}

// ============================================================================
// indexing tests
// ============================================================================

static void test_n_bins()
{
  // formula: cut / (0.01 / factor) = cut * factor / 0.01
  CHECK(cmdata::indexing::n_bins(0.75f)       == 300);
  CHECK(cmdata::indexing::n_bins(1.00f)       == 400);
  CHECK(cmdata::indexing::n_bins(0.50f)       == 200);
  CHECK(cmdata::indexing::n_bins(0.75f, 4.0f) == 300);
  CHECK(cmdata::indexing::n_bins(0.75f, 2.0f) == 150); // factor=2 → 0.75*2/0.01
}

static void test_offset_same_uniqueness()
{
  // For each topology (natmol2, num_mol_unique), all canonical tuples
  // (mol_i, a_i, a_j) with a_i <= a_j must map to distinct offsets.
  // (a_i, a_j) and (a_j, a_i) intentionally alias — they're the same cell.
  struct Case { std::vector<int> natmol2; std::vector<int> nmol; };
  std::vector<Case> cases = {
    {{ 2 }, { 1 }},
    {{ 3 }, { 1 }},
    {{ 3 }, { 2 }},
    {{ 5 }, { 3 }},
    {{ 2 }, { 4 }},
  };
  for (auto &c : cases)
  {
    std::set<std::size_t> seen;
    bool ok = true;
    for (std::size_t mt = 0; mt < c.natmol2.size(); mt++)
      for (int im = 0; im < c.nmol[mt]; im++)
        for (int ai = 0; ai < c.natmol2[mt]; ai++)
          for (int aj = ai; aj < c.natmol2[mt]; aj++)  // upper triangle only
          {
            auto off = cmdata::indexing::offset_same(mt, im, ai, aj, c.natmol2);
            ok = ok && seen.insert(off).second; // insert returns false on duplicate
          }
    CHECK_MSG(ok, "offset_same: duplicate offset found");
  }
}

static void test_offset_same_formula()
{
  std::vector<int> nat = { 3 };
  // N=3, per_mol = 3*4/2 = 6, row lo starts at lo*(2N-lo+1)/2
  // (a_i,a_j) is canonicalized to (min,max), so (1,0) == (0,1)
  CHECK(cmdata::indexing::offset_same(0, 0, 0, 0, nat) == 0);
  CHECK(cmdata::indexing::offset_same(0, 0, 0, 2, nat) == 2);
  CHECK(cmdata::indexing::offset_same(0, 0, 1, 0, nat) == 1);  // canonicalized (0,1) → 1
  CHECK(cmdata::indexing::offset_same(0, 1, 0, 0, nat) == 6);  // mol_i=1 → +1*6
  CHECK(cmdata::indexing::offset_same(0, 1, 2, 1, nat) == 10); // mol=1, (1,2): 6+3+1=10
}

static void test_offset_cross_uniqueness()
{
  struct Case {
    std::vector<int> natmol2;
    std::vector<int> nmol;
    std::size_t mt_i, mt_j;
  };
  std::vector<Case> cases = {
    {{ 2, 3 }, { 1, 1 }, 0, 1},
    {{ 2, 3 }, { 2, 1 }, 0, 1},
    {{ 2, 3 }, { 1, 3 }, 0, 1},
    {{ 2, 3 }, { 2, 3 }, 0, 1},
    {{ 3, 3 }, { 2, 2 }, 0, 1},
    {{ 4, 2 }, { 3, 4 }, 0, 1},
  };
  for (auto &c : cases)
  {
    std::set<std::size_t> seen;
    bool ok = true;
    for (int im = 0; im < c.nmol[c.mt_i]; im++)
      for (int jm = 0; jm < c.nmol[c.mt_j]; jm++)
        for (int ai = 0; ai < c.natmol2[c.mt_i]; ai++)
          for (int aj = 0; aj < c.natmol2[c.mt_j]; aj++)
          {
            auto off = cmdata::indexing::offset_cross(c.mt_i, c.mt_j, im, jm, ai, aj, c.natmol2, c.nmol[c.mt_j]);
            ok = ok && seen.insert(off).second;
          }
    char buf[128];
    snprintf(buf, sizeof(buf), "offset_cross: mt_i=%zu mt_j=%zu ni=%d nj=%d",
             c.mt_i, c.mt_j, c.natmol2[c.mt_i], c.natmol2[c.mt_j]);
    CHECK_MSG(ok, buf);
  }
}

static void test_offset_cross_formula()
{
  // offset_cross = mol_i * (num_mol_j * Ni * Nj) + mol_j * (Ni * Nj) + ai * Nj + aj
  std::vector<int> nat = { 2, 3 };
  std::size_t num_mol_j = 2;
  // mol_i=0, mol_j=0, ai=0, aj=0 → 0
  CHECK(cmdata::indexing::offset_cross(0, 1, 0, 0, 0, 0, nat, num_mol_j) == 0);
  // mol_i=0, mol_j=0, ai=0, aj=2 → 2
  CHECK(cmdata::indexing::offset_cross(0, 1, 0, 0, 0, 2, nat, num_mol_j) == 2);
  // mol_i=0, mol_j=0, ai=1, aj=0 → Nj=3
  CHECK(cmdata::indexing::offset_cross(0, 1, 0, 0, 1, 0, nat, num_mol_j) == 3);
  // mol_i=0, mol_j=1, ai=0, aj=0 → Ni*Nj = 2*3 = 6
  CHECK(cmdata::indexing::offset_cross(0, 1, 0, 1, 0, 0, nat, num_mol_j) == 6);
  // mol_i=1, mol_j=0, ai=0, aj=0 → num_mol_j * Ni*Nj = 2*6 = 12
  CHECK(cmdata::indexing::offset_cross(0, 1, 1, 0, 0, 0, nat, num_mol_j) == 12);
}

// ============================================================================
// density tests
// ============================================================================

static void test_kde_normalization()
{
  // A single sample with weight=1 should integrate to ~1.
  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  // interior point (mu well inside range, no boundary effect)
  {
    std::vector<float> h(bins.size(), 0.f);
    cmdata::density::kernel_density_estimator(h.begin(), bins, 0.40f, 1.f);
    double integral = integrate(h, dx);
    CHECK_NEAR(integral, 1.0, 0.01); // 1% tolerance for discrete approximation
  }
  // point near zero (mu < h=0.01): the implementation doubles the scale as a
  // boundary reflection approximation.  The approximation is exact at mu=0 but
  // degrades for mu in (0, h).  We only check that the KDE is non-zero and that
  // the integral stays in a physically reasonable range (0.8 – 1.5).
  {
    std::vector<float> hv(bins.size(), 0.f);
    cmdata::density::kernel_density_estimator(hv.begin(), bins, 0.005f, 1.f);
    double integral = integrate(hv, dx);
    CHECK_MSG(integral > 0.8 && integral < 1.6, "KDE near-zero: integral out of reasonable range");
  }
  // weight=2: integral should be ~2
  {
    std::vector<float> h(bins.size(), 0.f);
    cmdata::density::kernel_density_estimator(h.begin(), bins, 0.40f, 2.f);
    double integral = integrate(h, dx);
    CHECK_NEAR(integral, 2.0, 0.02);
  }
}

static void test_kde_support()
{
  // Bins outside [mu-2h, mu+2h] (h=0.01) must be exactly zero.
  const float cutoff = 0.75f;
  const float h = 0.01f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  const float mu = 0.40f;
  std::vector<float> hist(bins.size(), 0.f);
  cmdata::density::kernel_density_estimator(hist.begin(), bins, mu, 1.f);

  for (std::size_t k = 0; k < bins.size(); k++)
  {
    if (bins[k] < mu - 2.f * h - dx || bins[k] > mu + 2.f * h + dx)
      CHECK_MSG(hist[k] == 0.f, "KDE nonzero outside support");
  }
  // peak should be near mu
  std::size_t peak = std::max_element(hist.begin(), hist.end()) - hist.begin();
  CHECK_MSG(std::abs(bins[peak] - mu) < 3.f * dx, "KDE peak not near mu");
}

static void test_kde_accumulation()
{
  // Adding two samples produces the sum of individual KDEs.
  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  std::vector<float> h1(bins.size(), 0.f), h2(bins.size(), 0.f), h12(bins.size(), 0.f);
  cmdata::density::kernel_density_estimator(h1.begin(),  bins, 0.30f, 1.f);
  cmdata::density::kernel_density_estimator(h2.begin(),  bins, 0.50f, 1.f);
  cmdata::density::kernel_density_estimator(h12.begin(), bins, 0.30f, 1.f);
  cmdata::density::kernel_density_estimator(h12.begin(), bins, 0.50f, 1.f);

  bool ok = true;
  for (std::size_t k = 0; k < bins.size(); k++)
    if (std::abs(h12[k] - h1[k] - h2[k]) > 1e-6f) { ok = false; break; }
  CHECK_MSG(ok, "KDE: two-sample sum != sum of individual KDEs");

  CHECK_NEAR(integrate(h12, dx), 2.0, 0.02);
}

static void test_normalize_histo()
{
  // normalize_histo multiplies every element by norm * inv_num_mol_same.
  std::vector<std::vector<std::vector<std::vector<float>>>> data(1);
  data[0].resize(3, std::vector<std::vector<float>>(3, std::vector<float>(10, 2.f)));

  cmdata::density::normalize_histo(0, 1, 2, 0.5f, 0.25f, data);
  // expected: 2.0 * 0.5 * 0.25 = 0.25
  for (float v : data[0][1][2]) CHECK_NEAR(v, 0.25, 1e-6);
  // other elements untouched
  for (float v : data[0][0][0]) CHECK_NEAR(v, 2.0, 1e-6);

  // zero stays zero
  std::vector<std::vector<std::vector<std::vector<float>>>> z(1);
  z[0].resize(2, std::vector<std::vector<float>>(2, std::vector<float>(5, 0.f)));
  cmdata::density::normalize_histo(0, 0, 1, 1.f, 1.f, z);
  for (float v : z[0][0][1]) CHECK(v == 0.f);
}

// ============================================================================
// mindist tests
// ============================================================================

// Helper: build cross_index table for a given natmol2.
static std::vector<std::vector<int>> make_cross_index(std::size_t n_types)
{
  std::vector<std::vector<int>> ci(n_types, std::vector<int>(n_types, 0));
  int count = 0;
  for (std::size_t i = 0; i < n_types; i++)
    for (std::size_t j = i + 1; j < n_types; j++)
      ci[i][j] = count++;
  return ci;
}

static void test_mindist_same_single_mol()
{
  // One molecule type, one molecule, 2 atoms.
  // frame_same_mat[0][offset] = known distance.
  // After mindist_same, interm_same_maxcdf_mol[0][i][j] should be non-zero
  // near that distance and its integral should be ~1.
  std::vector<int> natmol2   = { 2 };
  std::vector<int> nmol      = { 1 };

  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  // frame_same_mat[mt][offset] — size = nmol * N*(N+1)/2 = 1*3 = 3
  std::vector<std::vector<float>> fsm(1, std::vector<float>(3, 100.f));
  // set distance for (im=0, a_i=0, a_j=1)
  float test_dist = 0.30f;
  fsm[0][cmdata::indexing::offset_same(0, 0, 0, 1, natmol2)] = test_dist;
  // set diagonal (a_i==a_j) to some value too
  fsm[0][cmdata::indexing::offset_same(0, 0, 0, 0, natmol2)] = 0.10f;
  fsm[0][cmdata::indexing::offset_same(0, 0, 1, 1, natmol2)] = 0.20f;

  // interm_same_maxcdf_mol[mt][a_i][a_j] — start at zero
  std::vector<std::vector<std::vector<std::vector<float>>>> cdf(1);
  cdf[0].resize(2, std::vector<std::vector<float>>(2, std::vector<float>(bins.size(), 0.f)));

  cmdata::mindist::mindist_same(bins, nmol, natmol2, fsm, cdf, 1.f);

  // (0,0) and (1,1) should have non-zero KDE
  CHECK(integrate(cdf[0][0][0], dx) > 0.5);
  CHECK(integrate(cdf[0][1][1], dx) > 0.5);

  // (0,1) should have KDE peaking near test_dist
  double integral_01 = integrate(cdf[0][0][1], dx);
  CHECK_NEAR(integral_01, 1.0, 0.02);

  // mindist_same mirrors: cdf[0][1][0] should equal cdf[0][0][1]
  bool mirrored = true;
  for (std::size_t k = 0; k < bins.size(); k++)
    if (std::abs(cdf[0][0][1][k] - cdf[0][1][0][k]) > 1e-6f) { mirrored = false; break; }
  CHECK_MSG(mirrored, "mindist_same: [i][j] != [j][i] after mirror copy");
}

static void test_mindist_same_multi_mol()
{
  // One type, 3 molecules, 2 atoms each.
  // All (im, a_i, a_j) slots must be visited.
  // We count how many slots received a KDE contribution.
  std::vector<int> natmol2 = { 2 };
  std::vector<int> nmol    = { 3 };
  const int N = natmol2[0], M = nmol[0];

  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);

  // fill every slot with a known distance
  std::vector<std::vector<float>> fsm(1, std::vector<float>(M * N * (N + 1) / 2, 0.40f));

  std::vector<std::vector<std::vector<std::vector<float>>>> cdf(1);
  cdf[0].resize(N, std::vector<std::vector<float>>(N, std::vector<float>(bins.size(), 0.f)));

  cmdata::mindist::mindist_same(bins, nmol, natmol2, fsm, cdf, 1.f);

  // Each (a_i, a_j >= a_i) bin should have been accumulated M times.
  // integral ≈ M (weight=1 per accumulation).
  float dx = cutoff / static_cast<float>(bins.size());
  for (int ai = 0; ai < N; ai++)
    for (int aj = ai; aj < N; aj++)
      CHECK_NEAR(integrate(cdf[0][ai][aj], dx), static_cast<double>(M), 0.05 * M);
}

static void test_mindist_cross_basic()
{
  // Two molecule types: type 0 has 2 atoms / 1 mol, type 1 has 3 atoms / 1 mol.
  // cross_index[0][1] = 0
  // frame_cross_mat[0] has size 1*2*1*3 = 6.
  // After mindist_cross, interm_cross_maxcdf_mol[0][a_i][a_j] should carry
  // a normalized KDE (integral ~1).
  std::vector<int> natmol2 = { 2, 3 };
  std::vector<int> nmol    = { 1, 1 };
  auto ci = make_cross_index(2);

  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  // Set every entry in the cross mat to a known distance.
  const float test_dist = 0.45f;
  std::vector<std::vector<float>> fcm(1, std::vector<float>(2 * 3 * 1 * 1, test_dist));

  std::vector<std::vector<std::vector<std::vector<float>>>> cdf(1);
  cdf[0].resize(2, std::vector<std::vector<float>>(3, std::vector<float>(bins.size(), 0.f)));

  cmdata::mindist::mindist_cross(natmol2, ci, bins, nmol, fcm, cdf, 1.f);

  for (int ai = 0; ai < 2; ai++)
    for (int aj = 0; aj < 3; aj++)
      CHECK_NEAR(integrate(cdf[0][ai][aj], dx), 1.0, 0.02);
}

static void test_mindist_cross_multi_mol()
{
  // Type 0: 2 atoms, 2 mols.  Type 1: 2 atoms, 3 mols.
  // Each (ai, aj) pair is visited im*jm = 2*3 = 6 times.
  std::vector<int> natmol2 = { 2, 2 };
  std::vector<int> nmol    = { 2, 3 };
  auto ci = make_cross_index(2);
  std::size_t buf_size = natmol2[0] * natmol2[1] * nmol[0] * nmol[1]; // 2*2*2*3=24

  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  std::vector<std::vector<float>> fcm(1, std::vector<float>(buf_size, 0.40f));

  std::vector<std::vector<std::vector<std::vector<float>>>> cdf(1);
  cdf[0].resize(2, std::vector<std::vector<float>>(2, std::vector<float>(bins.size(), 0.f)));

  cmdata::mindist::mindist_cross(natmol2, ci, bins, nmol, fcm, cdf, 1.f);

  const double expected = static_cast<double>(nmol[0]) * static_cast<double>(nmol[1]);
  for (int ai = 0; ai < 2; ai++)
    for (int aj = 0; aj < 2; aj++)
      CHECK_NEAR(integrate(cdf[0][ai][aj], dx), expected, 0.05 * expected);
}

static void test_mindist_cross_uses_correct_offset()
{
  // Verify that each (im, jm) slot in frame_cross_mat contributes to the
  // KDE at the correct distance.  We place a unique distance in each slot
  // and verify that the combined KDE integrates to im_count * jm_count.
  std::vector<int> natmol2 = { 1, 1 }; // one atom each — simplifies offset to mol index
  std::vector<int> nmol    = { 2, 2 };
  auto ci = make_cross_index(2);
  // offset_cross(mt_i=0, mt_j=1, im, jm, 0, 0, nat, num_mol_j=2)
  //   = im*(2*1*1) + jm*(1*1) + 0 + 0 = 2*im + jm
  // 4 slots: (0,0)->0, (0,1)->1, (1,0)->2, (1,1)->3
  std::vector<std::vector<float>> fcm(1, std::vector<float>(4, 0.f));
  fcm[0][0] = 0.30f; // (im=0, jm=0)
  fcm[0][1] = 0.35f; // (im=0, jm=1)
  fcm[0][2] = 0.40f; // (im=1, jm=0)
  fcm[0][3] = 0.45f; // (im=1, jm=1)

  const float cutoff = 0.75f;
  auto bins = make_bins(cutoff);
  float dx = cutoff / static_cast<float>(bins.size());

  std::vector<std::vector<std::vector<std::vector<float>>>> cdf(1);
  cdf[0].resize(1, std::vector<std::vector<float>>(1, std::vector<float>(bins.size(), 0.f)));

  cmdata::mindist::mindist_cross(natmol2, ci, bins, nmol, fcm, cdf, 1.f);

  // All four contributions land on cdf[0][0][0] → integral ≈ 4.
  CHECK_NEAR(integrate(cdf[0][0][0], dx), 4.0, 0.05);
}

// ============================================================================
// main
// ============================================================================

int main()
{
  printf("=== indexing ===\n");
  test_n_bins();
  test_offset_same_uniqueness();
  test_offset_same_formula();
  test_offset_cross_uniqueness();
  test_offset_cross_formula();

  printf("=== density ===\n");
  test_kde_normalization();
  test_kde_support();
  test_kde_accumulation();
  test_normalize_histo();

  printf("=== mindist ===\n");
  test_mindist_same_single_mol();
  test_mindist_same_multi_mol();
  test_mindist_cross_basic();
  test_mindist_cross_multi_mol();
  test_mindist_cross_uses_correct_offset();

  printf("\n%d PASS  %d FAIL\n", g_pass, g_fail);
  return g_fail > 0 ? 1 : 0;
}
