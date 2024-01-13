#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
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
        ~Simulate() {}
        void set_address_bits(int address_bits) { this->address_bits = address_bits;}
        void set_block_size(int block_size) {
            this->block_size = block_size;
            offset_bits = log2(block_size);
            cache_bits = address_bits - offset_bits;
        }
        void set_cache_sets(int cache_sets) {
            this->cache_sets = cache_sets;
            index_bits_num = log2(cache_sets);
            index_bits.reserve(index_bits_num);
        }
        void set_associativity(int associativity) { this->associativity = associativity;}
        int get_address_bits() { return address_bits;}
        int get_block_size() { return block_size;}
        int get_cache_sets() { return cache_sets;}
        int get_associativity() { return associativity;}
        void simulation(ofstream& outfile);
        void initialize(ifstream& infile, ofstream& outfile, int LSB_en = 0);
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
        vector<int> index_bits;
        string benchmark;
        map<string,Cache> simulate_cache; 
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
    simulate.initialize(infile, outfile, 0);
    simulate.simulation(outfile);
    infile.close();
    return 0;
}

void Simulate::initialize(ifstream& infile, ofstream& outfile, int LSB_en) {
    outfile << "Offset bit count: " << offset_bits << endl;
    outfile << "Indexing bit count: " << index_bits_num << endl;
    getline(infile, benchmark);
    string line;
    while(getline(infile, line)) {
        if(line == ".end") break;
        reverse(line.begin(), line.end());
        ref.push_back(line);
    }
    if(LSB_en){
        for(int i = 0; i < index_bits_num; i++) {
            index_bits.push_back(i);
        }
    }
    else{
        vector< vector<double> > correlation_matrix (cache_bits, vector<double>(cache_bits, 0));
        for(int i = 0; i < ref.size(); i++) {
            for(int j = 0; j < cache_bits; j++) {
                for(int k = j + 1; k < cache_bits; k++) {
                    if(ref[i][j+offset_bits] == ref[i][k+offset_bits]) {
                        correlation_matrix[j][k]++;
                        correlation_matrix[k][j]++;
                    }
                }
            }
        }
        for(int i = 0; i < cache_bits; i++) {
            for(int j = 0; j < cache_bits; j++) {
                if(correlation_matrix[i][j] > ref.size() - correlation_matrix[i][j]) 
                    correlation_matrix[i][j] = (double) (ref.size() - correlation_matrix[i][j]) / correlation_matrix[i][j];
                else
                    correlation_matrix[i][j] = (double) correlation_matrix[i][j] / (ref.size() - correlation_matrix[i][j]);
            }
        }
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
        for(int i=0; i<index_bits_num; i++) {
            int max_index = 0;
            for(int j = 0; j < cache_bits; j++) {
                if(quality[j] > quality[max_index]) max_index = j;
            }
            index_bits.push_back(max_index);
            for(int j = 0; j < cache_bits; j++) {
                quality[j] *= correlation_matrix[max_index][j];
            }
        }
        sort(index_bits.begin(), index_bits.end());
        index_bits.erase(unique(index_bits.begin(), index_bits.end()), index_bits.end());
        while(index_bits.size() < index_bits_num) {
            int flag = 0;
            for(int i = 0; i < index_bits.size();i++){
                if(index_bits[i] != i) {
                    index_bits.insert(index_bits.begin() + i, i);
                    flag = 1;
                    break;
                }
            }
            if(flag == 0) index_bits.push_back(index_bits.size());
        }
    }
    outfile << "Indexing bits:";
    for(int i = index_bits_num-1; i >= 0; i--) {
        outfile << " " << offset_bits + index_bits[i];
    }
    outfile << endl << endl;
}

void Simulate::simulation(ofstream& outfile) {
    int total_miss = 0;
    outfile << benchmark << endl;
    for(int i=0; i<ref.size(); i++){
        string address = ref[i].substr(offset_bits, address_bits - offset_bits);
        string index = "";
        string tag = "";
        for(int j = 0; j < cache_bits; j++) {
            if(find(index_bits.begin(), index_bits.end(), j) != index_bits.end()) {
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
            total_miss++;
            outfile << ref[i] <<  " " << "miss" << endl;
        }
    }
    outfile << ".end" << endl << endl;
    outfile << "Total cache miss count: " << total_miss << endl;
}