// Unit tests for cmdata::partition and cmdata::indexing::offset_cross.
//
// No GROMACS, no external test framework — plain C++17 with assertions and
// stderr reporting. Exit code 0 = all tests pass.
//
// Build standalone (no cmake required):
//   g++ -std=c++17 -Wall -I../src test_partition.cpp -o test_partition && ./test_partition

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "cmdata/indexing.hpp"
#include "cmdata/partition.hpp"

// ─────────────────────────────────────────────────────────────────────────────
// Minimal test harness
// ─────────────────────────────────────────────────────────────────────────────

static int g_pass = 0, g_fail = 0;

#define CHECK(expr)                                                            \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::fprintf(stderr, "FAIL  %s:%d  %s\n", __FILE__, __LINE__, #expr);   \
      ++g_fail;                                                                \
    } else {                                                                   \
      ++g_pass;                                                                \
    }                                                                          \
  } while (0)

#define CHECK_MSG(expr, msg)                                                   \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::fprintf(stderr, "FAIL  %s:%d  %s  (%s)\n",                         \
                   __FILE__, __LINE__, #expr, (msg).c_str());                  \
      ++g_fail;                                                                \
    } else {                                                                   \
      ++g_pass;                                                                \
    }                                                                          \
  } while (0)

// ─────────────────────────────────────────────────────────────────────────────
// Tuple helpers
// ─────────────────────────────────────────────────────────────────────────────

using Op6 = std::tuple<int,int,int,int,int,int>; // mt_i, mt_j, im, jm, i, j
using Op4 = std::tuple<int,int,int,int>;          // mt_i, im,   i,  j

// ─────────────────────────────────────────────────────────────────────────────
// Enumerate every cross op in canonical order  (mt_i, mt_j>mt_i, im, jm, i, j)
// This is the order mindist_cross iterates.
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<Op6> all_cross_ops(
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique
)
{
  std::vector<Op6> ops;
  for (int mt_i = 0; mt_i < (int)natmol2.size(); mt_i++)
    for (int mt_j = mt_i + 1; mt_j < (int)natmol2.size(); mt_j++)
      for (int im = 0; im < num_mol_unique[mt_i]; im++)
        for (int jm = 0; jm < num_mol_unique[mt_j]; jm++)
          for (int i = 0; i < natmol2[mt_i]; i++)
            for (int j = 0; j < natmol2[mt_j]; j++)
              ops.emplace_back(mt_i, mt_j, im, jm, i, j);
  return ops;
}

// ─────────────────────────────────────────────────────────────────────────────
// Enumerate every same op in canonical order  (mt_i, im, i, j>=i)
// This is the order mindist_same iterates.
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<Op4> all_same_ops(
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique
)
{
  std::vector<Op4> ops;
  for (int mt_i = 0; mt_i < (int)natmol2.size(); mt_i++)
    for (int im = 0; im < num_mol_unique[mt_i]; im++)
      for (int i = 0; i < natmol2[mt_i]; i++)
        for (int j = i; j < natmol2[mt_i]; j++)
          ops.emplace_back(mt_i, im, i, j);
  return ops;
}

// ─────────────────────────────────────────────────────────────────────────────
// Simulate mindist_cross for a single thread — returns the ops it would process
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<Op6> simulate_cross(
  const cmdata::indexing::CrossThreadIndices &idx,
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique
)
{
  std::vector<Op6> ops;
  if (idx.n_loop_operations_cross == 0) return ops;

  bool first_mtj = true, first_im = true, first_jm = true;
  bool first_i = true, first_j = true;
  int counter = 0;

  for (int mt_i = (int)idx.start_mti_cross; mt_i < (int)natmol2.size(); mt_i++)
  {
    for (int mt_j = first_mtj ? (int)idx.start_mtj_cross : mt_i + 1;
         mt_j < (int)natmol2.size(); mt_j++)
    {
      for (int im = first_im ? (int)idx.start_im_cross : 0;
           im < num_mol_unique[mt_i]; im++)
      {
        for (int jm = first_jm ? (int)idx.start_jm_cross : 0;
             jm < num_mol_unique[mt_j]; jm++)
        {
          for (int i = first_i ? (int)idx.start_i_cross : 0;
               i < natmol2[mt_i]; i++)
          {
            for (int j = first_j ? (int)idx.start_j_cross : 0;
                 j < natmol2[mt_j]; j++)
            {
              ops.emplace_back(mt_i, mt_j, im, jm, i, j);
              ++counter;
              if (counter == idx.n_loop_operations_cross) return ops;
              first_j = false;
            }
            first_j = false;
          }
          first_i = false;
        }
        first_jm = false;
      }
      first_im = false;
    }
    first_mtj = false;
  }
  return ops;
}

// ─────────────────────────────────────────────────────────────────────────────
// Simulate mindist_same for a single thread
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<Op4> simulate_same(
  const cmdata::indexing::SameThreadIndices &idx,
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique
)
{
  std::vector<Op4> ops;
  if (idx.n_loop_operations_same == 0) return ops;

  bool first_im = true, first_i = true, first_j = true;
  long int counter = 0;

  for (int mt_i = (int)idx.start_mti_same; mt_i < (int)natmol2.size(); mt_i++)
  {
    for (int im = first_im ? (int)idx.start_im_same : 0;
         im < num_mol_unique[mt_i]; im++)
    {
      for (int i = first_i ? (int)idx.start_i_same : 0; i < natmol2[mt_i]; i++)
      {
        for (int j = first_j ? (int)idx.start_j_same : i; j < natmol2[mt_i]; j++)
        {
          ops.emplace_back(mt_i, im, i, j);
          ++counter;
          if (counter == idx.n_loop_operations_same) return ops;
          first_j = false;
        }
        first_j = false;
      }
      first_i = false;
    }
    first_im = false;
  }
  return ops;
}

// ─────────────────────────────────────────────────────────────────────────────
// Core checker: partition coverage and no-overlap for cross case
// ─────────────────────────────────────────────────────────────────────────────
static void check_cross_partition(
  const std::string &label,
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique,
  int num_threads
)
{
  if (natmol2.size() < 2) return; // no cross case possible

  auto result = cmdata::partition::compute_thread_indices(
    natmol2, num_mol_unique, num_threads, /*do_same=*/false, /*do_cross=*/true
  );

  const auto expected = all_cross_ops(natmol2, num_mol_unique);
  const int total_ops = (int)expected.size();

  // Build a map from op → thread index that claimed it, to detect duplicates
  std::map<Op6, int> op_to_thread;
  std::vector<Op6> all_claimed;

  for (int tid = 0; tid < num_threads; tid++)
  {
    auto ops = simulate_cross(result.cross[tid], natmol2, num_mol_unique);

    // n_loop_operations == claimed count
    CHECK_MSG(
      (int)ops.size() == result.cross[tid].n_loop_operations_cross,
      label + " tid=" + std::to_string(tid) + " ops.size() vs n_loop"
    );

    for (auto &op : ops)
    {
      auto [it, inserted] = op_to_thread.emplace(op, tid);
      CHECK_MSG(inserted,
        label + " duplicate op (" +
        std::to_string(std::get<0>(op)) + "," +
        std::to_string(std::get<1>(op)) + "," +
        std::to_string(std::get<2>(op)) + "," +
        std::to_string(std::get<3>(op)) + "," +
        std::to_string(std::get<4>(op)) + "," +
        std::to_string(std::get<5>(op)) + ") claimed by tid=" +
        std::to_string(it->second) + " AND tid=" + std::to_string(tid)
      );
      all_claimed.push_back(op);
    }
  }

  // Total claimed == expected total
  CHECK_MSG(
    (int)all_claimed.size() == total_ops,
    label + " total ops: got " + std::to_string(all_claimed.size()) +
    " expected " + std::to_string(total_ops)
  );

  // Claimed set == expected set (i.e. no gaps)
  std::set<Op6> claimed_set(all_claimed.begin(), all_claimed.end());
  std::set<Op6> expected_set(expected.begin(), expected.end());
  CHECK_MSG(claimed_set == expected_set, label + " claimed set != expected set");

  // Verify the ops are delivered in the same ORDER as the canonical enumeration
  // when all thread op-lists are concatenated in tid order.
  std::vector<Op6> concat;
  for (int tid = 0; tid < num_threads; tid++)
  {
    auto ops = simulate_cross(result.cross[tid], natmol2, num_mol_unique);
    concat.insert(concat.end(), ops.begin(), ops.end());
  }
  CHECK_MSG(concat == expected,
    label + " concatenated order does not match canonical order"
  );
}

// ─────────────────────────────────────────────────────────────────────────────
// Core checker: partition coverage for same case
// ─────────────────────────────────────────────────────────────────────────────
static void check_same_partition(
  const std::string &label,
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique,
  int num_threads
)
{
  // same mode only makes sense when at least one type has >1 molecule
  bool any_multi = false;
  for (auto n : num_mol_unique) if (n > 1) { any_multi = true; break; }
  if (!any_multi) return;

  auto result = cmdata::partition::compute_thread_indices(
    natmol2, num_mol_unique, num_threads, /*do_same=*/true, /*do_cross=*/false
  );

  const auto expected = all_same_ops(natmol2, num_mol_unique);
  const int total_ops = (int)expected.size();

  std::map<Op4, int> op_to_thread;
  std::vector<Op4> all_claimed;

  for (int tid = 0; tid < num_threads; tid++)
  {
    auto ops = simulate_same(result.same[tid], natmol2, num_mol_unique);
    CHECK_MSG(
      (int)ops.size() == result.same[tid].n_loop_operations_same,
      label + " same tid=" + std::to_string(tid) + " ops.size() vs n_loop"
    );
    for (auto &op : ops)
    {
      auto [it, inserted] = op_to_thread.emplace(op, tid);
      CHECK_MSG(inserted, label + " same: duplicate op claimed by tid=" +
        std::to_string(it->second) + " AND tid=" + std::to_string(tid));
      all_claimed.push_back(op);
    }
  }

  CHECK_MSG(
    (int)all_claimed.size() == total_ops,
    label + " same total ops: got " + std::to_string(all_claimed.size()) +
    " expected " + std::to_string(total_ops)
  );

  std::vector<Op4> concat;
  for (int tid = 0; tid < num_threads; tid++)
  {
    auto ops = simulate_same(result.same[tid], natmol2, num_mol_unique);
    concat.insert(concat.end(), ops.begin(), ops.end());
  }
  CHECK_MSG(concat == expected, label + " same order mismatch");
}

// ─────────────────────────────────────────────────────────────────────────────
// Check offset_cross uniqueness: every (mol_i, mol_j, a_i, a_j) must map to a
// distinct slot in [0, num_mol_unique[i]*num_mol_unique[j]*natmol2[i]*natmol2[j])
// ─────────────────────────────────────────────────────────────────────────────
static void check_offset_cross(
  const std::string &label,
  const std::vector<int> &natmol2,
  const std::vector<int> &num_mol_unique
)
{
  for (int mt_i = 0; mt_i < (int)natmol2.size(); mt_i++)
  {
    for (int mt_j = mt_i + 1; mt_j < (int)natmol2.size(); mt_j++)
    {
      const std::size_t capacity =
        static_cast<std::size_t>(num_mol_unique[mt_i])
        * static_cast<std::size_t>(num_mol_unique[mt_j])
        * static_cast<std::size_t>(natmol2[mt_i])
        * static_cast<std::size_t>(natmol2[mt_j]);

      std::set<std::size_t> seen_offsets;
      bool all_unique = true, all_in_bounds = true;

      for (int im = 0; im < num_mol_unique[mt_i]; im++)
      {
        for (int jm = 0; jm < num_mol_unique[mt_j]; jm++)
        {
          for (int a_i = 0; a_i < natmol2[mt_i]; a_i++)
          {
            for (int a_j = 0; a_j < natmol2[mt_j]; a_j++)
            {
              std::size_t off = cmdata::indexing::offset_cross(
                mt_i, mt_j, im, jm, a_i, a_j, natmol2
              );
              if (off >= capacity) all_in_bounds = false;
              if (!seen_offsets.insert(off).second) all_unique = false;
            }
          }
        }
      }

      std::string pair = label + " cross(" + std::to_string(mt_i) + "," +
                         std::to_string(mt_j) + ") nmu=[" +
                         std::to_string(num_mol_unique[mt_i]) + "," +
                         std::to_string(num_mol_unique[mt_j]) + "]";

      CHECK_MSG(all_in_bounds, pair + ": offset out of bounds");
      CHECK_MSG(all_unique,    pair + ": aliased offsets (bug in formula)");
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Test topologies
// ─────────────────────────────────────────────────────────────────────────────

struct Topo {
  std::string label;
  std::vector<int> natmol2;
  std::vector<int> num_mol_unique;
};

static const std::vector<Topo> TOPOS = {
  // trivial: 1 type, single molecule — no cross, no same
  {"1t_1mol",       {5},        {1}},

  // 1 type, multi-molecule: same only
  {"1t_3mol",       {4},        {3}},
  {"1t_1mol_big",   {172},      {1}},
  {"1t_2mol_big",   {172},      {2}},

  // 2 types, 1 molecule each: minimal cross
  {"2t_1each",      {3,4},      {1,1}},

  // 2 types, multi-molecule: cross + same
  {"2t_2each",      {3,4},      {2,2}},
  {"2t_3_2",        {5,3},      {3,2}},

  // 3 types — the case that exposed the bug
  {"3t_1each",      {2,2,2},    {1,1,1}},
  {"3t_3_1_2",      {2,2,2},    {3,1,2}},  // concrete example from the bug trace
  {"3t_vary",       {4,3,2},    {2,3,1}},
  {"3t_protein_lg", {172,6,3},  {2,4,8}},  // protein + ligand + water (scaled)
  {"3t_1multi",     {3,3,3},    {1,3,1}},  // only middle type has multi-mol

  // 4 types: stress test
  {"4t_all1",       {2,3,4,2},  {1,1,1,1}},
  {"4t_mixed",      {4,2,3,2},  {2,1,3,2}},
};

static const std::vector<int> THREAD_COUNTS = {1, 2, 3, 4, 7, 8, 16};

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main()
{
  std::printf("Running partition unit tests...\n\n");

  // ── cross partition tests ──────────────────────────────────────────────────
  std::printf("=== Cross partition coverage tests ===\n");
  for (const auto &topo : TOPOS)
  {
    if (topo.natmol2.size() < 2) continue;
    for (int nt : THREAD_COUNTS)
    {
      std::string lbl = topo.label + "/threads=" + std::to_string(nt);
      check_cross_partition(lbl, topo.natmol2, topo.num_mol_unique, nt);
    }
  }

  // ── same partition tests ───────────────────────────────────────────────────
  std::printf("=== Same partition coverage tests ===\n");
  for (const auto &topo : TOPOS)
  {
    for (int nt : THREAD_COUNTS)
    {
      std::string lbl = topo.label + "/threads=" + std::to_string(nt);
      check_same_partition(lbl, topo.natmol2, topo.num_mol_unique, nt);
    }
  }

  // ── offset_cross uniqueness tests ─────────────────────────────────────────
  std::printf("=== offset_cross uniqueness tests ===\n");
  for (const auto &topo : TOPOS)
  {
    if (topo.natmol2.size() < 2) continue;
    check_offset_cross(topo.label, topo.natmol2, topo.num_mol_unique);
  }

  // ── summary ───────────────────────────────────────────────────────────────
  std::printf("\n");
  if (g_fail == 0)
    std::printf("ALL %d checks PASSED.\n", g_pass);
  else
    std::printf("%d FAILED / %d PASSED (total %d)\n", g_fail, g_pass, g_pass + g_fail);

  return g_fail > 0 ? 1 : 0;
}
