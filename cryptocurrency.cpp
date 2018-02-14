#include<iostream>
#include<vector>
#include<random>
#include<deque>
#include<map>
#include<math.h>
#include<time.h>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <fstream>

using namespace std;


/// Global Variables ///

int N = 20; //Number of peers
float txn_mean = 0.2; // Mean of exponential transaction distribution function
float cpu_power_mean = 150.0; // Mean time of a block generated in the system
float z = 0.4; // Probability of a fast node
float ff = 0.1; // Probability of fast-fast node connection
float fs = 0.5; // Probability of fast-slow node connection
float ss = 0.2; // Probability of slow-slow node connection

int txn_counter = 0; // purpose : maintain unique txn ids
int block_counter = 1; // purpose : maintain unique block ids

/// Randomizor Functions ///

float exp_dist(float mean) { // Exponential Distribution
	random_device rd;
	mt19937 gen(rd());
	exponential_distribution<> dist(1.0 / mean);
	return dist(gen);
}

float uni_dist(float start, float end) { // Uniform Distribution
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dist(start, end);
	return dist(gen);
}


/// Structures ///

struct tnx { // Transaction Structure
	int id;
  	int send_id;
  	int recv_id;
  	float amount;
  	double time;

  	tnx(int txn_id, int id1, int id2, float val, double time) { // Create a new transaction
  		id 		= txn_id;
  		send_id = id1;
  		recv_id = id2;
  		amount  = val;
  		this->time = time;
  	}
};

struct block { // Block Structure
	int peer_id;
	int prevblockID;
	int blockID;
	double time_arrival;
	vector<tnx> spent;
};

struct event { // Global Time Queue Event Element Structure
	int event_type = 0; // 0 for tranx, 1 for block generation, 2 for add block 
	int peer_id;
	double time;
	double time_gen = 0.0; 
	block * blk; 

	event(int event_type, int peer_id, double time) { // Create a new event
		this->event_type = event_type;
		this->peer_id = peer_id;
		this->time = time;
		blk = NULL;
	}
};

struct node { // Wrapper for Blockchain Block
	block *blk;
	node *parent;
	vector<node *> nextBlocks;
	vector<float> peer_amount;

	node() { // Create a new node
		peer_amount.assign(N, 100.0);
	}		
};	

struct ans_long_chain { // Wrapper for obtaining Head of longest chain
	node *last_node;
	double time;
	int length;

	// Constructors
	ans_long_chain(double time, int length) {
		this->time = time;
		this->length = length;
	}

	ans_long_chain(node *last_block, double time, int length) {
		this->last_node = last_block;
		this->time = time;
		this->length = length;
	}
};	

vector<event> time_simulator; // Global Time Queue for Simulator


/// Helper Functions ///

// Add element to Global Time Vector
void sorted_event_add(event ev) {
	if(time_simulator.size() == 0)
		time_simulator.push_back(ev);
	int i = 0;
	for(; i < time_simulator.size(); i++) {
		if(time_simulator[i].time > ev.time)
			break;
	}
	time_simulator.insert(time_simulator.begin()+i, ev);
}

// Add element to Peer Transaction Vector
void sorted_add(tnx trans, vector<tnx> &globalQueueTnx) {
	if(globalQueueTnx.size() == 0)
		globalQueueTnx.push_back(trans);
	int i = 0;
	for(; i < globalQueueTnx.size(); i++) {
		if(globalQueueTnx[i].time > trans.time)
			break;
	}
	globalQueueTnx.insert(globalQueueTnx.begin()+i, trans);
} 


/// Classes /// 

class blockchain { // The blockchain class
  public:
	node *root;

	blockchain() { // Create a new blockchain with Genisis block
		root = new node();
		root->blk = new block();
		root->blk->blockID = 0;
		root->parent = NULL;	
	}

	// Returns the node in the block chain with a certain block id
	node* find_block(node *start, int block_id) {
		if(start->blk->blockID == block_id)
			return start;

		for(auto child : start->nextBlocks) {
			node *ans = find_block(child, block_id);
			if(ans != NULL)
				return ans;
		}
		return NULL;
	}

	// Finds the lowest common ancestor between two leaf nodes
	node* find_common_anc(node *prev, node *curr) {
		// make list of all parents of previous node
		node *tmp = prev;
		unordered_set<node*> parents;
		while(tmp !=NULL) {
			parents.insert(tmp->parent);
			tmp = tmp->parent;
		}

		// find first parent of current node in list of parents
		node* tempo = curr;
		while(tempo) {
			if(parents.find(tempo) != parents.end())
				return tempo;
			tempo = tempo->parent; 
		}
		return NULL; // Should not happen
	}

	// Returns the head of the longest chain of the blockchain
	ans_long_chain getLongestChain(node *start, int length) {
		// If leaf node
		if(start->nextBlocks.size() == 0) {
			ans_long_chain ans(start->blk->time_arrival, length+1);
			ans.last_node = start;
			return ans;
		}

		// Scan all Children
		node* lastNode;
		int max_length = -1;
		double min_time = std::numeric_limits<double>::max();
		for(auto child : start->nextBlocks) {
			ans_long_chain child_ans = getLongestChain(child, length+1);
			if((child_ans.length > max_length) || ((child_ans.length == max_length) && (child_ans.time < min_time))) {
				max_length = child_ans.length;
				min_time = child_ans.time;
				lastNode = child_ans.last_node;
			}
		}
		ans_long_chain final_ans(lastNode, min_time, max_length);
		return final_ans;
	}

	// Roll Back on transactions in case of change of longest blockchain
	void update(vector <tnx> &globalQueueTnx, node* prev, node* curr) {
		node *anc = find_common_anc(prev, curr);

		// Add all the transactions to txn queue in the blocks being rolled back
		node *tmp = prev;
		while(tmp != anc) { 
			for(auto it : tmp->blk->spent)
				sorted_add(it, globalQueueTnx);
			tmp = tmp->parent;
		}

		// Remove all the transactions from txn queue in the blocks being added 
		tmp = curr;
		while(tmp != anc) {
			unordered_set<int> trans; 
			for(auto it : tmp->blk->spent)
				trans.insert(it.id);
			for(int i = 0; i < globalQueueTnx.size(); i++)
				if(trans.find(globalQueueTnx[i].id) != trans.end())
					globalQueueTnx.erase(globalQueueTnx.begin()+i);
			tmp = tmp->parent; 
		}
	}

	// Add new block to the blockchain
	node* add_block(int prevblockID, int blockID, double time_arrival) {
		node *lastNode = find_block(root, prevblockID); // Get head of the blockchain
		if(lastNode == NULL)
 			return NULL;
		node *newNode = new node();
		newNode->parent = lastNode;
		newNode->blk = new block();
		newNode->blk->prevblockID = prevblockID;
		newNode->blk->blockID = blockID;
		newNode->blk->time_arrival = time_arrival;
		lastNode->nextBlocks.push_back(newNode);
		return newNode;
	}

	node* generate_node(double time) {
		node *current = new node();
		current->blk = new block();
		current->blk->blockID = block_counter;
		ans_long_chain long_chain = getLongestChain(root, 0);
		node *longestChain = (long_chain.last_node);
		current->parent = longestChain;
		current->blk->prevblockID = longestChain->blk->blockID;
		current->blk->time_arrival = time;
		longestChain->nextBlocks.push_back(current);
		block_counter++;
		return current;
	}
};	

class network; // Global network class declaration

class peer { // The peer node class 
  private:
	int ID;
	int activation;
	int connected;
	vector <tnx> globalQueueTnx; // Unspent transactions
	vector<int> blocks_rec;
	unordered_map<int, vector<block>> orphan_blk;
	double lastBlockArrival = 0.0; 

  public:
  	bool type;
  	node *longest_one;
  	blockchain *chain;
  	float block_mean;

	peer(int id, bool speed, int active) { // Create a new peer for the network
		ID = id;
		type = speed;
		activation = active;
		connected = 0;
		chain = new blockchain();
		longest_one = chain->root;
		block_mean = exp_dist(cpu_power_mean);
	}

	int get_id() { // Get node id
		return ID;
	}

	bool get_type() { // Get Speed of node
		return type;
	}

	bool is_unactive() { // Check if node has been activated or not
		if (connected < activation) 
			return true;
		else
			return false;
	}

	void add_conn() {
		connected += 1;
	}

	bool txn_exists(int id) { // Check if a txn exists in the txn queue
		for(auto it : globalQueueTnx)
			if(it.id == id)
				return true;
		return false;	
	}

	bool blk_exist(int id) { // Check if block already received
		for(auto it : blocks_rec)
			if(it == id)
				return true;
		return false;	
	}

	void generate_transaction(double time, network * tmp); 

	void generate_block(double time, double time_gen, network * tmp);

	void add_txn(tnx trans, network * tmp); 
		
	void add_blk(block &blk, bool again);

	void add_blk_sim(block &blk, network * tmp);

};

// Maintains the cryptocurrency network
class network {
  private:
  	
	map<int, vector< pair< peer*, float > > > adjlist;
	deque<peer*> unactive;
	vector<int> visited;

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

  	vector<peer*> nodelist;
  	peer* findnode(int id) { // Find a node in the network
		int nsize = nodelist.size();
		for (int i = 0; i < nsize; i++) {
			if (nodelist[i]->get_id() == id)
				return nodelist[i];
		}
		return NULL; // return NULL if node does not exist
	}

	void addnode(int id, bool type, int active) { // Add a node to the network
		cout<<"Adding node : id = " << id << " type = " << type << " activation = " << active << endl;
		if (findnode(id) == NULL) {
			peer* newnode = new peer(id, type, active);
			nodelist.push_back(newnode);
			unactive.push_back(newnode);
		}
	}

	void addconn(peer* send_node, peer* recv_node) { // Add a connection between two nodes in network
		int send_id = send_node->get_id();
		int recv_id = recv_node->get_id();
		float latency = uni_dist(10, 500) / 1000;
		add_edge(send_id, recv_node, latency);
		add_edge(recv_id, send_node, latency);
		send_node->add_conn();
		recv_node->add_conn();
	}

	bool findconn(int u, int v) { // Find if connection exists between two nodes of the network
		int usize = adjlist[u].size();
		for (int i = 0; i < usize; i++) {
			if (adjlist[u][i].first->get_id() == v)
				return true;
		}
		return false;
	}

	float findlatency(int u, int v) { // Calculates the link latency between two nodes of the network 
		int usize = adjlist[u].size();
		for (int i = 0; i < usize; i++) {
			if (adjlist[u][i].first->get_id() == v)
				return adjlist[u][i].second;
		}
		return 0.0;
	}

	bool is_connected() { // Checks if all the nodes of the network are activated
		if (unactive.size() == 0)
			return true;
		else
			return false;
	}

	bool try_to_connect(bool type1, bool type2) { // Attempts to make a connection between two chosen nodes
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

	void connect_graph() { // Create a connected network of nodes
		// Add peers the network
		for (int i = 0; i < N; i++) {
			float node_type_prob = uni_dist(0,1);
			bool node_type = false;
			int activation;
			if (node_type_prob < z)
				node_type = true;
			if (node_type)
				activation = 1.5*log(N);
			else 
				activation = log(N); 
			addnode(i, node_type, activation);
		}

		// Connect Peers
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

	void dfs(int node_id) { // DFS helper function to see if network is strongly connected
		visited.push_back(node_id);
		for( const auto& peerlist : adjlist[node_id] ) {
			if(find(visited.begin(), visited.end(), peerlist.first->get_id())==visited.end())
				dfs(peerlist.first->get_id());
		}
	}

	bool disconnected() { // Check to see if the network is strongly connected
		dfs(0);
		if (visited.size() == N)
			return false;
		else
			return true;
	}

	void reset_graph() { // Creates a new network
		adjlist.clear();
		unactive.clear();
		visited.clear();
	}

	void print_graph() { // Helper function to print the graph and connections
		for( const auto& peerlist : adjlist ) {
			int node_id = peerlist.first;
			cout << "\nnode " << node_id << " : ";
			for (int i = 0; i < peerlist.second.size(); i++) {
				peer* node = peerlist.second[i].first;
				cout << node->get_id() << " ";
			}
		}
	}

	float get_latency(peer* send_node, peer*recv_node, int m_size) { // Calculate the total latency between two nodes in the network
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
		float total_latency = (float)m_size / c + exp_dist(d_mean) + p;
		return total_latency;
	}

	void broadcast_tnx(tnx &trans, peer* recv_node, int send_id) { // Broadcasts a txn to all other peers
		if (recv_node->txn_exists(trans.id))
			return;
		else {
			recv_node->add_txn(trans, this);
			for( const auto& peerlist : adjlist[recv_node->get_id()] ) {
				if(peerlist.first->get_id() != send_id)
					broadcast_tnx(trans, peerlist.first, recv_node->get_id());
			}
		}
	}

	void broadcast_blk(block &blk, peer* recv_node, int send_id) { // Broadcasts a block to all other peers
		if (recv_node->blk_exist(blk.blockID) )
			return;
		else {
			if(recv_node->get_id() != blk.peer_id)
				recv_node->add_blk_sim(blk, this);
			for( const auto& peerlist : adjlist[recv_node->get_id()] ) {
				if(peerlist.first->get_id() != send_id)
					broadcast_blk(blk, peerlist.first, recv_node->get_id());
			}
		}
	}	
	
};

// Generate a transaction at a peer
void peer::generate_transaction(double time, network * tmp) { // for intial transaction, time would be zero 
	float pay_amt = 0.2 * longest_one->peer_amount[this->ID] * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
	int rec = rand() % N;
	while(rec == ID)
		rec = rand() % N;
	tnx trans(++txn_counter, ID, rec, pay_amt, time);
	double next_time = time + exp_dist(txn_mean);
	event nextTrans(0, this->ID, next_time);
	sorted_event_add(nextTrans);
	tmp->broadcast_tnx(trans, this, -1);
}

// Generate a block at a peer
void peer::generate_block(double time, double time_gen, network * tmp) { 
	// add only those transaction for which time is less than this
	if(lastBlockArrival <= time && lastBlockArrival > time_gen)
		return;
	node *current = chain->generate_node(time);
	current->blk->peer_id = this->ID;
	for(int i = 0; i < globalQueueTnx.size(); i++) {
		if(globalQueueTnx[i].time <= time)
			current->blk->spent.push_back(globalQueueTnx[i]);
		else {
			globalQueueTnx.erase(globalQueueTnx.begin(), globalQueueTnx.begin()+i);
			break;
		}
	}
	tmp->broadcast_blk(*(current->blk), this, -1);
}

// Add a transaction to the peer txn list (adjusted for network latencies)
void peer::add_txn(tnx trans, network *tmp) { 
	if(trans.send_id != this->ID) {
		double latency = tmp->get_latency(tmp->nodelist[trans.send_id], this, 0);
		trans.time += latency;
	}	
	sorted_add(trans, globalQueueTnx);
}

void peer::add_blk_sim(block &blk, network *tmp) {
	blocks_rec.push_back(blk.blockID);
	if(blk.peer_id != this->ID) {
		double latency = tmp->get_latency(tmp->nodelist[blk.peer_id], this, 1); 
		blk.time_arrival += latency;
	}
	event blk_gen(2, this->ID, blk.time_arrival);
	blk_gen.blk = new block();
	blk_gen.blk->peer_id = blk.peer_id;
	blk_gen.blk->prevblockID = blk.prevblockID;
	blk_gen.blk->blockID = blk.blockID;
	blk_gen.blk->time_arrival = blk.time_arrival;
	blk_gen.blk->peer_id = blk.peer_id;
	for(auto it : blk.spent)
		blk_gen.blk->spent.push_back(it);
	sorted_event_add(blk_gen);
}

// Attempt to add a block to a peer blockchain
void peer::add_blk(block &blk, bool again) {
	node *current = chain->add_block(blk.prevblockID, blk.blockID, blk.time_arrival);
	bool flag = false;
	if(!current) {
		orphan_blk[blk.prevblockID].push_back(blk);
		flag = true;
	}	
	if(blk.peer_id != this->ID) {
		if(!again) {
			lastBlockArrival = blk.time_arrival;
			// generate a random variable and create event for next block gen
			event blk_gen(1, this->ID, lastBlockArrival+exp_dist(block_mean));
			blk_gen.time_gen = lastBlockArrival;
			sorted_event_add(blk_gen);
		}
		if(flag)
			return;
		unordered_set<int> trans_id;
		for(auto it : blk.spent)
			trans_id.insert(it.id);
		for(int i = 0; i < globalQueueTnx.size(); i++) {
			if(trans_id.find(globalQueueTnx[i].id) != trans_id.end()) {
				globalQueueTnx.erase(globalQueueTnx.begin() + i);
				i--;
			}
		}
	}
	if(flag)
		return;
	ans_long_chain long_chain = chain->getLongestChain(chain->root, 0);
	node* new_longest_one = long_chain.last_node;
	if(blk.prevblockID != longest_one->blk->blockID && new_longest_one->blk->blockID == blk.blockID)	
		chain->update(globalQueueTnx, longest_one, new_longest_one);
	longest_one = new_longest_one;
	for(int i = 0; i < N; i++)
		current->peer_amount[i] = current->parent->peer_amount[i];
	for(auto it : blk.spent) {
		current->blk->spent.push_back(it);
		if(it.amount > current->peer_amount[it.send_id])
			continue;
		current->peer_amount[it.send_id] -= it.amount;
		current->peer_amount[it.recv_id] += it.amount;
		current->peer_amount[current->blk->peer_id] += 50; // mining fees
	}
	if(orphan_blk.find(current->blk->blockID) != orphan_blk.end()) {
		for(auto it: orphan_blk[current->blk->blockID]) {
			this->add_blk(it, true);
		}
		orphan_blk.erase(current->blk->blockID);
	}
}

/// Visualizer ///

// Create edges of graph via graph traversal
void chain_span(node* root) {
	for (auto & child : root->nextBlocks) {
		cout << root->blk->blockID << " -- " << child->blk->blockID << endl;
		chain_span(child);
	}
}

// Function to create a graph file to visualize local blockchain of a peer
void create_visual(string file, blockchain* object) {
	ofstream out(file);
	streambuf *coutbuf = cout.rdbuf();
	cout.rdbuf(out.rdbuf());
	cout << "graph {" << endl;
	chain_span(object->root);
	cout << "}" << endl;
	cout.rdbuf(coutbuf);
}

/// Main Function ///

int main() {

	// Create P2P network for cryptocurrency
	network mycoin;

	// Connect peers in the network
	mycoin.connect_graph();
	while(mycoin.disconnected()){
		mycoin.reset_graph();
		mycoin.connect_graph();
	}

	// Initialising the queue
	for(int i = 0; i < N; i++) {
		cout<<mycoin.nodelist[i]->block_mean<<", ";
		event txn(0, i, exp_dist(txn_mean));
		sorted_event_add(txn);
		event blk(1, i, exp_dist(mycoin.nodelist[i]->block_mean));
		sorted_event_add(blk);
	}
	cout<<endl;

	// Select certain nodes for visualization
	int fast_low_id = 0;
	double fast_low_pow = -1.0; 
	int fast_high_id = 0;
	double fast_high_pow = 100000.0;
	int slow_low_id = 0;
	double slow_low_pow = -1.0;
	int slow_high_id = 0;
	double slow_high_pow = 1000000.0;

	for(int i = 0; i < N; i++) {
		if(mycoin.nodelist[i]->type) {
			if(mycoin.nodelist[i]->block_mean >fast_low_pow) {
				fast_low_pow = mycoin.nodelist[i]->block_mean;
				fast_low_id = i;
			}
			if(mycoin.nodelist[i]->block_mean < fast_high_pow) {
				fast_high_pow = mycoin.nodelist[i]->block_mean;
				fast_high_id = i;
			}
		}
		else {
			if(mycoin.nodelist[i]->block_mean > slow_low_pow) {
				slow_low_pow = mycoin.nodelist[i]->block_mean;
				slow_low_id = i;
			}
			if(mycoin.nodelist[i]->block_mean < slow_high_pow) {
				slow_high_pow = mycoin.nodelist[i]->block_mean;
				slow_high_id = i;
			}
		}
		
	}

	int counter = 0 ;
	while(true) { // Run the simulation (Get events from the event queue)
		counter++;
		if(time_simulator.size() != 0) {	
			event current = time_simulator[0];
			time_simulator.erase(time_simulator.begin());
			if(current.event_type == 0)
				mycoin.nodelist[current.peer_id]->generate_transaction(current.time, &mycoin);
			else if(current.event_type ==1)
				mycoin.nodelist[current.peer_id]->generate_block(current.time, current.time_gen, &mycoin);	
  			else if(current.event_type ==2)
				mycoin.nodelist[current.peer_id]->add_blk(*(current.blk), false);	
  		
	  		if (current.time > 75.0) {
	  			cout <<" COUNTER " <<endl;
	  			create_visual("testing_20_f_l_p150.txt", mycoin.findnode(fast_low_id)->chain);
	  			create_visual("testing_20_s_l_p150.txt", mycoin.findnode(slow_low_id)->chain);
	  			create_visual("testing_20_f_h_p150.txt", mycoin.findnode(fast_high_id)->chain);
	  			create_visual("testing_20_s_h_p150.txt", mycoin.findnode(slow_high_id)->chain);
	  			for(int i = 0; i < N; i++)
					cout<<mycoin.nodelist[0]->longest_one->peer_amount[i]<<", ";
				cout<<endl;
	  			break;
	  		}
  	    }
	}
	return 0;
}


