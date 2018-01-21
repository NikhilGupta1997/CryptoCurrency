#include<iostream>
#include<vector>
#include<random>
#include<deque>
#include<map>
#include<math.h>

using namespace std;

// Global Variables
int N = 10; //Number of peers
float txn_mean = 100.0; // Mean of exponential transaction distribution function
float z = 0.5; // Probability of a fast node
float ff = 0.5;
float fs = 0.25;
float ss = 0.1;

// Randomizor Functions
float exp_dist(float mean) {
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<> dist(1.0 / mean);
	return dist(gen);
}

float uni_dist(float start, float end) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dist(start, end);
	return dist(gen);
}
	
// Classes
class peer {
  private:
	int ID;
	float amount;
	bool type;
	int activation;
	int connected;
	
  public:

	peer(int id, bool speed, int active) {
		ID = id;
		amount = 0;
		type = speed;
		activation = active;
		connected = 0;
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

	bool is_unactive() {
		if (connected < activation) 
			return true;
		else
			return false;
	}

	void add_conn() {
		connected += 1;
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
	map<int, vector< pair< peer*, float > > > adjlist;
	deque<peer*> unactive;

	void add_edge(int send_id, peer* recv_node, float latency) {
		if(adjlist.find(send_id) == adjlist.end()) {
			vector< pair< peer*, float > > vec;
			vec.push_back(make_pair(recv_node, latency));
			adjlist[send_id] = vec;
		}
		else 
			adjlist[send_id].push_back(make_pair(recv_node, latency));
	}

  public:

  	peer* findnode(int id) {
		int nsize = nodelist.size();
		for (int i = 0; i < nsize; i++) {
			if (nodelist[i]->get_id() == id)
				return nodelist[i];
		}
		return NULL; ///// return NULL if node does not exist
	}

	void addnode(int id, bool type, int active) {
		cout<<"Adding node : id = " << id << " type = " << type << " activation = " << active << endl;
		if (findnode(id) == NULL) {
			peer* newnode = new peer(id, type, active);
			nodelist.push_back(newnode);
			unactive.push_back(newnode);
		}
	}

	void addconn(peer* send_node, peer* recv_node) {
		int send_id = send_node->get_id();
		int recv_id = recv_node->get_id();
		float latency = uni_dist(10, 500) / 1000;
		add_edge(send_id, recv_node, latency);
		add_edge(recv_id, send_node, latency);
		send_node->add_conn();
		recv_node->add_conn();
	}

	bool findconn(int u, int v) {
		int usize = adjlist[u].size();
		for (int i = 0; i < usize; i++) {
			if (adjlist[u][i].first->get_id() == v)
				return true;
		}
		return false;
	}

	float findlatency(int u, int v) {
		int usize = adjlist[u].size();
		for (int i = 0; i < usize; i++) {
			if (adjlist[u][i].first->get_id() == v)
				return adjlist[u][i].second;
		}
		return 0.0;
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

	bool is_connected() {
		if (unactive.size() == 0)
			return true;
		else
			return false;
	}

	bool try_to_connect(bool type1, bool type2) {
		float conn_prob;
		if (type1 != type2) 
			conn_prob = fs;
		else if (type1) 
			conn_prob = ff;
		else
			conn_prob = ss;
		float prob = uni_dist(0,1);
		if (prob < conn_prob)
			return true;
		else
			return false;
	}

	void connect_graph() {
		while(is_connected() == false) {
			peer* node = unactive.front();
			bool node_type = node->get_type();
			int node_id = node->get_id();
			while(node->is_unactive()) {

				peer* conn_node;
				bool connected;
				do {
					int node_num = rand() % nodelist.size();
					conn_node = nodelist[node_num];
					if (findconn(node_id, conn_node->get_id())) {
						continue;
					}
					connected = try_to_connect(node_type, conn_node->get_type());
				} while(conn_node->get_id() == node_id || connected == false || findconn(node_id, conn_node->get_id()));
				addconn(node, conn_node);
			}
			unactive.pop_front();
		}
	}

	void print_graph() {
		for( const auto& peerlist : adjlist ) {
			int node_id = peerlist.first;
			cout << "\nnode " << node_id << " : ";
			for (int i = 0; i < peerlist.second.size(); i++) {
				peer* node = peerlist.second[i].first;
				cout << node->get_id() << " ";
			}
		}
	}

	float get_latency(peer* send_node, peer*recv_node, int m_size) {
		bool type1 = send_node->get_type();
		bool type2 = recv_node->get_type();
		int c;
		if (type1 != type2) 
			c = 5;
		else if (type1) 
			c = 100;
		else
			c = 5;
		float d_mean = 96.0 * 1000 / (c * 1024 * 1024);
		float p = findlatency(send_node->get_id(), recv_node->get_id());
		cout << "\np = " << p << endl;
		cout << "c = " << c << endl;
		cout << "d_mean = " << d_mean << endl;
		float total_latency = (float)m_size / c + exp_dist(d_mean) + p;
		return total_latency;
	}
};

int main() {
	// Create P2P network for cryptocurrency
	network mycoin;
	
	// Add peers the network
	for (int i = 0; i < N; i++) {
		float node_type_prob = uni_dist(0,1);
		bool node_type = false;
		int activation;
		if (node_type_prob < z)
			node_type = true;
		if (node_type)
			activation = sqrt(N);
		else 
			activation = log(N); 
		mycoin.addnode(i, node_type, activation);
	}

	// Connect peers in the network
	mycoin.connect_graph();

	mycoin.print_graph();

	cout << mycoin.get_latency(mycoin.findnode(0), mycoin.findnode(2), 1);

	cout<<"\nend"<<endl;
	return 0;
}