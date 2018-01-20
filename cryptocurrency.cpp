#include<iostream>
#include<vector>
#include<random>

using namespace std;

// Global Variables
int N = 10; //Number of peers
float txn_mean = 100.0; // Mean of exponential transaction distribution function
float z = 0.5; // Probability of a fast node
	
// Classes
class peer {
  private:
	int ID;
	float amount;
	bool type;
	
  public:
	vector<int> conn;

	peer(int id, bool speed) {
		ID = id;
		amount = 0;
		type = speed;
	}

	int get_id() { // Get node id
		return ID;
	}

	float get_amount() {
		return amount;
	}

	void set_amount(float val) {
		amount = val;
	}

	void add_amount(float val) {
		amount += val;
	}

	void sub_amount(float val) { 
		amount += val;
	}

	bool get_type() { // Get Speed of node
		return type;
	}

	void add_conn(int id) { // Add peer to list of connected nodes
		conn.push_back(id);
	}
};

struct tnx {
  public:
  	int send_id;
  	int recv_id;
  	float amount;

  	tnx(int id1, int id2, float val) { // Create a new transaction
  		send_id = id1;
  		recv_id = id2;
  		amount = val;
  	}
};

// Maintains the cryptocurrency network
class network {
  private:
  	vector<peer*> nodelist;
	vector<peer*> *adjlist;

  public:
  	network() {
  		adjlist = new vector<peer*>[N];
  	}

  	peer* findnode(int id) {
		int nsize = nodelist.size();
		for (int i = 0; i < nsize; i++) {
			if (nodelist[i]->get_id() == id)
				return nodelist[i];
		}
		return NULL; ///// return NULL if node does not exist
	}

	void addnode(int id, bool type) {
		if (findnode(id) == NULL) {
			peer* newnode = new peer(id, type);
			nodelist.push_back(newnode);
		}
	}

	void addconn(int u, int v) {
		adjlist[u].push_back(findnode(v));
	}

	bool findconn(int u, int v) {
		int usize = adjlist[u].size();
		for (int i = 0; i < usize; i++) {
			if (adjlist[u][i]->get_id() == v)
				return true;
		}
		return false;
	}

	bool validate(tnx transaction) {
		int send_id = transaction.send_id;
		int txn_amount = transaction.amount;
		peer* send_node = findnode(send_id);
		if (send_node->get_amount() < txn_amount) 
			return false;
		else
			return true;
	}
};


int main() {
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<> exp_dist(10);
	uniform_real_distribution<> uni_dist(0, 1);
	return 0;
}