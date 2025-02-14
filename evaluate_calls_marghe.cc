#include <bits/stdc++.h>
#include <err.h>
#include <sysexits.h>
#include <math.h>

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

 /* while (file) {
    file >> name >> true_taxon;
    if(true_taxon > 0) {
    	truth_map[name] = true_taxon;
    }
  }*/

   while (getline(file, line)) {
	istringstream iss(line);
   	 iss >> name >> true_taxon;
	 if(true_taxon > 0) {
		truth_map[name] = true_taxon;
	 }
  }

}

void UpdateMapElements(int tax_id, unordered_map<int, int> &ump) {

  if (ump.find(tax_id) != ump.end()) {
      ump[tax_id] += 1; // increment map's value for key 'tax_id'
    } else { // key not found
      ump[tax_id] = 1;
    }
}

void CountTruth(unordered_map<int, int> &truth_count, unordered_map<string, int> &truth_map, const unordered_map<int, int> &parent_map, 
  const unordered_map<int, string> &rank_map, const string rank) {

  int correct = 0;

  for (auto el : truth_map) {
    correct = el.second;

    while (correct > 0) {
      if (rank_map.find(correct) == rank_map.end()) {
      	correct = 0;
      	break;
      }
      if (rank_map.at(correct) == rank) {
        break;
      }
      correct = parent_map.at(correct);
    }

    if(correct > 0) {
	   UpdateMapElements(correct, truth_count);
    }
  }
}
int not_in_tree = 0;

// is a an ancestor of b?
bool isAncestor(const unordered_map<int, int> &parent_map, int a, int b) {
  while (b > 0) {
    if (a == b)
      return true;
    if (parent_map.find(b) != parent_map.end()) {
      b = parent_map.at(b);
    } else {
      not_in_tree++;
      return false;
    }
  }
  return false;
}

double Mean(const unordered_map<int, int> &ump, const size_t n) {
  int sum = 0;

  for (auto el : ump) {
    sum += el.second;
  }

  return (1.0 * sum) / n;
}

double PearsonCorrelation(const unordered_map<int,int> &class_count, const unordered_map<int, int> &truth_count, const size_t n) {

 /* size_t n = truth_count.size();
  double truth_mean = Mean(truth_count, n);
  double class_mean = Mean(class_count, n);
  double num_sum = 0;
  double sum_square_truth = 0;
  double sum_square_class = 0;

  for (auto el : class_count) {
      num_sum += el.second * truth_count.at(el.first);
      sum_square_class += pow(1.0 * el.second, 2.0);
      sum_square_truth += pow(truth_count.at(1.0 * el.first), 2.0);
  }

  return (1.0 * num_sum - n * truth_mean * class_mean) / (sqrt(sum_square_class - n * pow(class_mean, 2.0)) * sqrt(sum_square_truth - n * pow(truth_mean, 2.0)));*/

  double sum_prod = 0;
  double class_sum = 0;
  double truth_sum = 0;
  double sum_square_truth = 0;
  double sum_square_class = 0;

  for (auto el : class_count) {
  	sum_prod += 1.0 * el.second * truth_count.at(el.first);
  	class_sum += 1.0 * el.second;
  	sum_square_class += 1.0 * el.second * el.second;
  	truth_sum += 1.0 * truth_count.at(el.first);
  	sum_square_truth += 1.0 * truth_count.at(el.first) * truth_count.at(el.first);
  }

/* cout << n << endl;
  cout << sum_prod <<endl;
  cout << class_sum << endl;
  cout << truth_sum << endl;
  cout << sum_square_class << endl;
  cout << sum_square_truth << endl;*/

  double numerator = n * sum_prod - class_sum  * truth_sum;
  double radicand_1 = n * sum_square_class - class_sum * class_sum;
  double radicand_2 = n * sum_square_truth - truth_sum * truth_sum;
  double radicand = radicand_1 * radicand_2;
  double denominator = sqrt(radicand);

/*  cout << radicand_1 << endl;
  cout << radicand_2 << endl;
  cout << numerator << endl;
  cout << denominator << endl;*/

  return numerator / denominator;

}

void ScoreCalls(const string &calls_file, const string &rank, const unordered_map<int, string> &rank_map,
  const unordered_map<int, int> &parent_map, unordered_map<string, int> &truth_map, unordered_map<int, int> &truth_count)
{
  ifstream file(calls_file);
  if (! file)
    err(EX_NOINPUT, "error opening %s", calls_file.c_str());

  int tp = 0; //true positive
  int fp = 0; //false positive
  int fn = 0; //false negative
  int ok = 0;  // correct but above rank
  int no = 0;  // no defined true value at rank
  set<int> tp_set; // true positive set

  unordered_map<int, int> class_count;
  size_t n = truth_map.size();
  string code, name;
  int taxon;
  string line;
  while (getline(file, line)) {
    istringstream iss(line);
    iss >> name >> taxon;
    if (truth_map.count(name) == 0) { // skip calls for things not in truth set
      continue;
    }
    if (taxon == 0) {  // unclassified, make FN
      truth_map.erase(name);
      fn++;
    }
    else {
      int correct = truth_map.at(name);
      truth_map.erase(name);
      if (rank_map.count(correct) == 0)  // handle cases where correct taxon isn't there
        correct = 0;
      // elevate correct taxon to specified rank
      while (correct > 0) {
        if (rank_map.at(correct) == rank) {
          break;
        }
        correct = parent_map.at(correct);
      }
      if (! correct) {
        no++;  // ran off tree without finding rank!
      }
      else {
        if (isAncestor(parent_map, correct, taxon)) {
          UpdateMapElements(correct, class_count);
          tp++;
            tp_set.insert(correct);

        }
        else {
          if (isAncestor(parent_map, taxon, correct))
            ok++;
          else {
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
  double pearson = PearsonCorrelation(class_count, truth_count, n);
  //printf("Analysis at %s level:\n", rank.c_str());
  printf("TruePositive\tFalsePositive\tFalseNegative\tCorrectAboveRank\tNoDefinedTrueValueAtRank\tTotalCount\tSensitivity\tPrecision\tF1\tPearson\tDistinctTPAtRank\tNotInTree\n");
  printf("%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%zu\t%d\n", tp, fp, fn, ok, no, (tp+fp+fn+ok+no), sens, prec, f1, pearson, tp_set.size(), not_in_tree);
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

  unordered_map<int,int> truth_count;
  CountTruth(truth_count, truth_map,  parent_map, rank_map, rank);

  ScoreCalls(calls_file, rank, rank_map, parent_map, truth_map, truth_count);

  return 0;
}
