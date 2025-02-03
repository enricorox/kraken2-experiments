#include <bits/stdc++.h>
#include <err.h>
#include <sysexits.h>

using namespace std;

void PopulateTaxonomy(unordered_map<int, string> &rank_map, unordered_map<int, int> &parent_map, const string &nodes_file) {
  ifstream file(nodes_file);
  if (! file)
    err(EX_NOINPUT, "error opening %s", nodes_file.c_str());
  string line, rank;
  int node_id = 0, parent_id = 0;

  const string delim = "\t|\t";
  while (getline(file, line)) {
    line.pop_back();
    line.pop_back();
    size_t pos1, pos2;
    pos1 = 0;
    int field_ct = 0;
    bool finished = false;
    while (field_ct++ < 10 && ! finished) {
      pos2 = line.find(delim, pos1);
      string token;
      if (pos2 == string::npos) {
        token = line.substr(pos1);
        finished = true;
      }
      else {
        token = line.substr(pos1, pos2 - pos1);
        pos1 = pos2 + delim.size();
      }

      switch (field_ct) {
        case 1:
          node_id = (int) stoul(token);
          break;
        case 2:
          parent_id = (int) stoul(token);
          break;
        case 3:
          rank = token;
          finished = true;
          break;
      }
    }  // end tokenizing loop
    if (node_id == 1)
      parent_id = 0;
    parent_map[node_id] = parent_id;
    rank_map[node_id] = rank;
  }
}

void ReadTruth(unordered_map<string, int> &truth_map, const string &truth_file) {
  ifstream file(truth_file);
  if (! file)
    err(EX_NOINPUT, "error opening %s", truth_file.c_str());
  string line, name;
  int true_taxon;

  while (file) {
    file >> name >> true_taxon;
    truth_map[name] = true_taxon;
  }
}

int not_in_tree = 0;

// is a_up an ancestor of b_down?
bool isAncestor(const unordered_map<int, int> &parent_map, int a_up, int b_down) {
  if(parent_map.find(b_down) == parent_map.end()) { // handle if not in tree
      not_in_tree++;
      return false;
  }

  while (b_down > 0) {
    if (a_up == b_down)
      return true;
    b_down = parent_map.at(b_down);
  }
  return false;
}

void ScoreCalls(const string &calls_file, const string &rank, const unordered_map<int, string> &rank_map,
  const unordered_map<int, int> &parent_map, unordered_map<string, int> &truth_map)
{
  ifstream file(calls_file);
  if (! file)
    err(EX_NOINPUT, "error opening %s", calls_file.c_str());

  int tp = 0;
  int fp = 0;
  int fn = 0;
  int ok = 0;  // correct but above rank
  int no = 0;  // no defined true value at rank
  set<int> tp_set;

  string code, name;
  int taxon;
  string line;

  while (getline(file, line)) {
    istringstream iss(line);
    iss >> name >> taxon;

    if (truth_map.count(name) == 0)  // skip calls for things not in truth set
      continue;

    if (taxon == 0) {  // unclassified, make FN
      truth_map.erase(name);
      fn++;
    }else { // classified in some way
      int correct = truth_map.at(name); // correct taxa for the read 'name'
      truth_map.erase(name);
      if (rank_map.count(correct) == 0)  // handle cases where correct taxon isn't there
        correct = 0;

      // elevate correct taxon to user specified rank
      while (correct > 0) {
        if (rank_map.at(correct) == rank)
          break;
        correct = parent_map.at(correct);
      }

      if (correct == 0) { // root reached
        no++;  // ran off tree without finding rank!
      }
      else { // internal node or leaf
        if (isAncestor(parent_map, correct, taxon)) { // correct
            tp++;
            tp_set.insert(taxon);
        }
        else {
          if (isAncestor(parent_map, taxon, correct)) // correct above rank
            ok++;
          else { // completely wrong branch
            fp++;
          }
        }
      }
    }
  }
  fn += truth_map.size();  // any uncalled frags go here
  double sens = (tp * 1.0) / (tp + fn + fp + ok);  // includes FP+OK because we need full set of eval'd frags in denominator
  double prec = (tp * 1.0) / (tp + fp);
  double f1 = sens * prec == 0 ? 0.0 : 2 * sens * prec / (sens + prec);
  printf("TruePositive\tFalsePositive\tFalseNegative\tCorrectAboveRank\tNoDefinedTrueValueAtRank\tTotalCount\tSensitivity\tPrecision\tF1\tDistinctTPAtRank\tNotInTree\n");
  printf("%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%zu\t%d\n", tp, fp, fn, ok, no, (tp+fp+fn+ok+no), sens, prec, f1, tp_set.size(), not_in_tree);
}

int main(int argc, char **argv) {
  if (argc != 5) {
    errx(EX_USAGE, "Usage: evaluate_calls <nodes.dmp> <rank> <truth.tsv> <calls.tsv>");
  }
  string nodes_file = argv[1];
  string rank = argv[2];
  string truth_file = argv[3];
  string calls_file = argv[4];

  unordered_map<int, string> rank_map;
  unordered_map<int, int> parent_map;
  PopulateTaxonomy(rank_map, parent_map, nodes_file);

  unordered_map<string, int> truth_map;
  ReadTruth(truth_map, truth_file);

  ScoreCalls(calls_file, rank, rank_map, parent_map, truth_map);

  return 0;
}
