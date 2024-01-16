#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <set>
using namespace std;

class Block{
    public:
    int NRU;
    string tag;
    Block() {
        NRU = 1;
        tag = "";
    }
    ~Block() {
    }
    bool hit(string tag) {
        if(this->tag == tag) {
            NRU = 0;
            return true;
        }
        return false;
    }
};

class Cache{
    public:
        Cache() {
            associativity = 0;
        }
        Cache(int associativity) {
            this->associativity = associativity;
            if(associativity > 0) blocks.resize(associativity);
        }
        ~Cache() {}
        bool hit(string tag) {
            for(int i = 0; i < associativity; i++) {
                if(blocks[i].hit(tag)) return true;
            }
            for(int i = 0; i < associativity; i++) {
                if(blocks[i].NRU == 1) {
                    blocks[i].tag = tag;
                    blocks[i].NRU = 0;
                    return false;
                }
            }
            for(int i = 0; i < associativity; i++) {
                blocks[i].NRU = 1;
            }
            blocks[0].tag = tag;
            blocks[0].NRU = 0;
            return false;
        }
    private:
        vector<Block> blocks;
        int associativity;
};


class Simulate{
    public:
        Simulate() {
            address_bits = 0;
            block_size = 0;
            cache_sets = 0;
            associativity = 0;
        }
        ~Simulate() {
            if(correlation_matrix != NULL) {
                delete correlation_matrix;
            } 
        }
        void set_address_bits(int address_bits) { this->address_bits = address_bits;}
        void set_block_size(int block_size) {
            this->block_size = block_size;
            offset_bits = log2(block_size);
            cache_bits = address_bits - offset_bits;
        }
        void set_cache_sets(int cache_sets) {
            this->cache_sets = cache_sets;
            index_bits_num = log2(cache_sets);
        }
        void set_associativity(int associativity) { this->associativity = associativity;}
        int get_address_bits() { return address_bits;}
        int get_block_size() { return block_size;}
        int get_cache_sets() { return cache_sets;}
        int get_associativity() { return associativity;}
        void simulation();
        void recursion(int bits_remain, vector<double> quality, set<int> index_bits);
        void initialize(ifstream& infile);
        void set_corr_matrix(ifstream& infile, ofstream& outfile);
        void output(ofstream& outfile);
    private:
        int address_bits;
        int block_size;
        int cache_sets;
        int associativity;
        int offset_bits;
        int index_bits_num;
        int cache_bits;
        vector<string> ref;
        vector<string> tags;
        set< set < int > > index_bits;
        string benchmark;
        map<string,Cache> simulate_cache; 
        vector< vector< double > >* correlation_matrix;
        int miss = 1e9;
        set < int > best = set < int >();
};

int main(int argc, char *argv[]) {
    ifstream infile;
    ofstream outfile;
    Simulate simulate;
    infile.open(argv[1]);
    if (!infile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    outfile.open(argv[3], ios::out);
    if (!outfile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    string line;
    for(int i = 0; i < 4; i++) {
        int tmp;
        infile >> line >> tmp;
        switch (i){
        case 0:
            simulate.set_address_bits(tmp);
            line.replace(line.find("_"), 1, " ");
            outfile << line << " " << simulate.get_address_bits() << endl;
            break;
        case 1:
            simulate.set_block_size(tmp);
            line.replace(line.find("_"), 1, " ");
            outfile << line << " " << simulate.get_block_size() << endl;
            break;
        case 2:
            simulate.set_cache_sets(tmp);
            line.replace(line.find("_"), 1, " ");
            outfile << line << " " << simulate.get_cache_sets() << endl;
            break;
        case 3:
            simulate.set_associativity(tmp);
            outfile << line << " " << simulate.get_associativity() << endl;
            break;
        }
    }
    outfile << endl;
    infile.close();
    infile.open(argv[2]);
    simulate.set_corr_matrix(infile, outfile);
    simulate.initialize(infile);
    simulate.simulation();
    simulate.output(outfile);
    infile.close();
    return 0;
}

void Simulate::initialize(ifstream& infile) {
    vector<double> quality (cache_bits, 0);
    for(int i = 0; i < ref.size(); i++) {
        for(int j = 0; j < cache_bits; j++) {
            quality[j] += (ref[i][j+offset_bits] == '1') ? 1 : 0;
        }
    }
    for(int i = 0; i < cache_bits; i++) {
        if(quality[i] > (ref.size() - quality[i])) 
            quality[i] = (double) (ref.size() - quality[i]) / quality[i];
        else
            quality[i] = (double) quality[i] / (ref.size() - quality[i]);
    }
    recursion(index_bits_num, quality, set<int>());
}

void Simulate::recursion(int bits_remain, vector<double> quality, set<int> index_bits){
    if(bits_remain == 0) {
        int i = 0;
        while (index_bits.size() < index_bits_num) {
            index_bits.insert(i++);
        }
        this->index_bits.insert(index_bits);
        return;
    }
    double max = 0;
    for(int j = 0; j < cache_bits; j++) {
        if(quality[j] > max) max = quality[j];
    }
    int flag = 0;
    for(vector<double>::iterator it = quality.begin(); it != quality.end(); ++it) {
        if(*it > 0.2){
            flag = 1;
            int max_index = it - quality.begin();
            vector<double> tmp = quality;
            for(int j = 0; j < cache_bits; j++) {
                tmp[j] *= (*correlation_matrix)[max_index][j];
            }
            set<int> tmp2 = index_bits;
            tmp2.insert(max_index);
            recursion(bits_remain - 1, tmp, tmp2);
        }
    }
    if(flag == 0) {
        for(vector<double>::iterator it = quality.begin(); it != quality.end(); ++it) {
            if(*it == max){
                int max_index = it - quality.begin();
                vector<double> tmp = quality;
                for(int j = 0; j < cache_bits; j++) {
                    tmp[j] *= (*correlation_matrix)[max_index][j];
                }
                set<int> tmp2 = index_bits;
                tmp2.insert(max_index);
                recursion(bits_remain - 1, tmp, tmp2);
            }
        }
    }
}

void Simulate::simulation() {
    for(set< set < int > >::iterator it = index_bits.begin(); it != index_bits.end(); ++it) {
        int tmp_miss = 0;
        for(int i=0; i<ref.size(); i++){
            string address = ref[i].substr(offset_bits, address_bits - offset_bits);
            string index = "";
            string tag = "";
            for(int j = 0; j < cache_bits; j++) {
                if((*it).find(j) != (*it).end()) {
                    index += address[j];
                }
                else {
                    tag += address[j];
                }
            }
            if(simulate_cache.find(index) == simulate_cache.end()) {
                simulate_cache[index] = Cache(associativity);
            }
            bool hit = simulate_cache[index].hit(tag);
            if(!hit) {
                tmp_miss++;
            }
        }
        if(tmp_miss < miss) {
            miss = tmp_miss;
            best = *it;
        }
        simulate_cache.clear();
    }
}

void Simulate::set_corr_matrix(ifstream& infile, ofstream& outfile) {
    outfile << "Offset bit count: " << offset_bits << endl;
    outfile << "Indexing bit count: " << index_bits_num << endl;
    getline(infile, benchmark);
    string line;
    while(getline(infile, line)) {
        if(line == ".end") break;
        reverse(line.begin(), line.end());
        ref.push_back(line);
    }
    correlation_matrix = new vector< vector<double> >(cache_bits, vector<double>(cache_bits, 0));
    for(int i = 0; i < ref.size(); i++) {
        for(int j = 0; j < cache_bits; j++) {
            for(int k = j + 1; k < cache_bits; k++) {
                if(ref[i][j+offset_bits] == ref[i][k+offset_bits]) {
                    (*correlation_matrix)[j][k]++;
                    (*correlation_matrix)[k][j]++;
                }
            }
        }
    }
    for(int i = 0; i < cache_bits; i++) {
        for(int j = 0; j < cache_bits; j++) {
            if((*correlation_matrix)[i][j] > ref.size() - (*correlation_matrix)[i][j]) 
                (*correlation_matrix)[i][j] = (double) (ref.size() - (*correlation_matrix)[i][j]) / (*correlation_matrix)[i][j];
            else
                (*correlation_matrix)[i][j] = (double) (*correlation_matrix)[i][j] / (ref.size() - (*correlation_matrix)[i][j]);
        }
    }
}

void Simulate::output(ofstream& outfile) {
    outfile << "Indexing bits:";
    for(set<int>::iterator i = best.begin(); i != best.end(); ++i) {
        outfile << " " << *i + offset_bits;
    }
    outfile << endl << endl;
    outfile << benchmark << endl;
    int total_miss = 0;
    for(int i=0; i<ref.size(); i++){
        string address = ref[i].substr(offset_bits, address_bits - offset_bits);
        string index = "";
        string tag = "";
        for(int j = 0; j < cache_bits; j++) {
            if(best.find(j) != best.end()) {
                index += address[j];
            }
            else {
                tag += address[j];
            }
        }
        if(simulate_cache.find(index) == simulate_cache.end()) {
            simulate_cache[index] = Cache(associativity);
        }
        bool hit = simulate_cache[index].hit(tag);
        reverse(ref[i].begin(), ref[i].end());
        if(hit) {
            outfile << ref[i] <<  " " << "hit" << endl;
        }
        else {
            outfile << ref[i] <<  " " << "miss" << endl;
        }
    }
    outfile << ".end" << endl << endl;
    outfile << "Total cache miss count: " << miss << endl;
}