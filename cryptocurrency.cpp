#include<iostream>
#include<vector>
#include<random>

using namespace std;

// Global Variables
int N = 10; //Number of peers
float txn_mean = 100.0; // Mean of exponential transaction distribution function

// Random Generation Machine
random_device rd;
mt19937 gen(rd());
	
// Classes
class peer {
  private:
	int ID;
	float amount;
	bool type;
	exponential_distribution<> dist(txn_mean);
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

	float get_txn_time() { // inter-arrival times between transactions
		return dist(gen)
	}
};

class tnx {
  private:
  	int send_id;
  	int recv_id;
  	float amount;

  	tnx(int id1, int id2, float val) { // Create a new transaction
  		send_id = id1;
  		recv_id = id2;
  		amount = val;
  	}
};

int main() {
	return 0;
}