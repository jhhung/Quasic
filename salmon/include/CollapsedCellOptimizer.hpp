#ifndef COLLAPSED_CELL_OPTIMIZER_HPP
#define COLLAPSED_CELL_OPTIMIZER_HPP

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <thread>

#include <boost/filesystem.hpp>

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"
#include "GZipWriter.hpp"
#include "EquivalenceClassBuilder.hpp"

#include "AlevinOpts.hpp"
#include "SingleCellProtocols.hpp"
#include "WhiteList.hpp"
#include "DedupUMI.hpp"
#include "TranscriptGroup.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"
#include "concurrentqueue.h"

#include <boost/math/special_functions/digamma.hpp>
#include <boost/random.hpp>

namespace bfs = boost::filesystem;
using JqueueT = moodycamel::ConcurrentQueue<uint32_t>;
using eqMapT = libcuckoo::cuckoohash_map<TranscriptGroup, SCTGValue, TranscriptGroupHasher>;
using tgrouplabelt = std::vector<uint32_t>;
using tgroupweightvec = std::vector<double>;
using SCExpT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;
using EqMapT = libcuckoo::cuckoohash_map<TranscriptGroup, SCTGValue, TranscriptGroupHasher>;

constexpr double digammaMin = 1e-10;

struct CellState {
  bool inActive;
};

struct CellTypePar {
  bool doCelltypeQuan = false;
  bool hasgeneBlacklist = false;
	double gamma = 0;
  double ARI = 0;
  double converge_ARI = 0.975;
  double pca_dim = 15; 
  double rare_gene_cut_off = 0.2;
	std::vector<bool> inGeneBlacklist;

  bool write_sig_file = false;
  bool write_iteration_count = false;
  bool writr_fileter_cells = true;
  bool writr_cluster_score = true;
  

	std::string write_iteration_file_name = "";
  std::string write_fileter_cells_file_name = "";
  std::string write_cluster_score_file_name = "";
};

struct Cell_cluster {
	uint32_t cluster;
	std::string name;
	std::vector<double> pca_distance;
	std::vector<double> score;
  double cluster_score;
  bool is_confidence;
  Cell_cluster(){};
	Cell_cluster(uint32_t clu, const std::string& s) : cluster(clu), name(s) {};

	bool operator > (const Cell_cluster& cell) const {
		return (name < cell.name);
	}
};

std::vector<Cell_cluster> read_cluster_file(std::string file_name);

class CollapsedCellOptimizer {

public:
  using VecType = std::vector<std::atomic<double>>;
  using SerialVecType = std::vector<double>;
  CollapsedCellOptimizer();

  template <class ProtocolT>
  bool optimize(EqMapT& freqMap,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                AlevinOpts<ProtocolT>& aopt,
                GZipWriter& gzw,
                std::vector<std::string>& trueBarcodes,
                std::vector<uint32_t>& umiCount,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
};

bool runPerCellEM(double& totalNumFrags, size_t numGenes,
                  CollapsedCellOptimizer::SerialVecType& alphas,
                  const CollapsedCellOptimizer::SerialVecType& priorAlphas,
                  std::vector<SalmonEqClass>& salmonEqclasses,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  bool initUniform, bool useVBEM, uint32_t cur_iter, 
                  const std::vector<double>& Cell_profile,
                  const std::vector<double>& Cell_type_profile, 
                  const std::vector<bool>& rare_gene_label,
                  CellTypePar& CTP, const Cell_cluster& cur_cell);

void optimizeCell(std::vector<std::string>& trueBarcodes,
                  const std::vector<std::vector<double>>& priorAlphas,
                  std::atomic<uint32_t>& barcode,
                  size_t totalCells, uint32_t umiEditDistance, eqMapT& eqMap,
                  std::deque<std::pair<TranscriptGroup, uint32_t>>& orderedTgroup,
                  std::shared_ptr<spdlog::logger>& jointlog,
                  std::vector<uint32_t>& umiCount,
                  std::vector<CellState>& skippedCB,
                  bool verbose, GZipWriter& gzw, bool noEM, bool useVBEM,
                  bool quiet, std::atomic<double>& totalDedupCounts,
                  std::atomic<uint32_t>& totalExpGeneCounts, double priorWeight,
                  spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                  uint32_t numGenes, uint32_t umiLength,
                  uint32_t numBootstraps, uint32_t numGibbsSamples,
                  bool naiveEqclass, bool dumpUmiGraph,
                  bool dumpCellEq, bool useAllBootstraps,
                  bool initUniform, CFreqMapT& freqCounter, bool dumpArboFragCounts,
                  spp::sparse_hash_set<uint32_t>& mRnaGenes,
                  spp::sparse_hash_set<uint32_t>& rRnaGenes,
                  std::atomic<uint64_t>& totalUniEdgesCounts,
                  std::atomic<uint64_t>& totalBiEdgesCounts, 
                  uint32_t cur_iter, std::map<std::string, std::vector<double>>& GEP, 
                  const std::vector<std::vector<double>>& Sig_mat, 
                  const std::vector<std::vector<double>>& Sig_mat_proportion, 
                  spp::sparse_hash_map<std::string, uint32_t>& cellIdxcluster, 
                  CellTypePar& CTP, std::vector<Cell_cluster>& y_cells);

using VecT = CollapsedCellOptimizer::SerialVecType;

#endif // COLLAPSED_CELL_OPTIMIZER_HPP
