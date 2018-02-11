#include<iostream>
#include<vector>
#include<random>
#include<deque>
#include<map>
#include<math.h>
#include<time.h>
#include <limits>
using namespace std;

// Global Variables
int N = 10; //Number of peers
float txn_mean = 100.0; // Mean of exponential transaction distribution function
float block_mean = 500.0; // TODO : Nikhil
float z = 0.5; // Probability of a fast node
float ff = 0.5;
float fs = 0.25;
float ss = 0.1;


int txn_counter = 0;
int block_counter = 1; // purpose : maintain unique block ids
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

struct event
{
	int event_type = 0; //0 for tranx, 1 for block
	int peer_id;
	double time;
	event(int event_type, int peer_id, double time)
	{
		this->event_type = event_type;
		this->peer_id = peer_id;
		this->time = time;
	}
};

vector<event> time_simulator;

void sorted_event_add(event ev)
{
	if(time_simulator.size() == 0)
		time_simulator.push_back(ev);
	int i = 0;
	for( ; i < time_simulator.size(); i++)
	{
		if(time_simulator[i].time > ev.time)
			break;
	}
	time_simulator.insert(time_simulator.begin()+i, ev);
}



struct tnx {
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
  		time = this->time;
  	}
};

struct block
{
	int peer_id;
	int prevblockID;
	int blockID;
	double time_arrival;
	vector<tnx> unspent;
	//TODO :  I should have a list of transactions as well or make a map for every block and transactions in it
};
struct node
{
	block *blk;
	node *parent;
	vector<node *> nextBlocks;
	vector<float> peer_amount;
	vector<tnx> unspent_trans;
	node()
	{
		peer_amount.assign(N, 100.0);
	}		
};	
struct ans_long_chain
{
	node *last_node;
	double time;
	int length;
	ans_long_chain(double time, int length)
	{
		this->time = time;
		this->length = length;
	}
	ans_long_chain(node *last_block, double time, int length)
	{
		this->last_node = last_block;
		this->time = time;
		this->length = length;
	}
};	
class blockchain
{
	
	node *root;

	node* find_block(node *start, int block_id)
	{
		if(start->blk->blockID == block_id)
		{
			return start;
		}
		// if(start->nextBlocks.size() == 0)
		// 	return NULL;
		for(auto child : start->nextBlocks)
		{
			node *ans = find_block(child, block_id);
			if(ans != NULL)
				return ans;
		}
		return NULL;
	}

	ans_long_chain getLongestChain(node *start, int length)
	{
		if(start->nextBlocks.size() == 0)
		{
			ans_long_chain ans(start->blk->time_arrival, length+1);
			ans.last_node = start;
			return ans;
		}
		node* lastNode;
		int max_length = -1;
		double min_time = std::numeric_limits<double>::max();
		for(auto child : start->nextBlocks)
		{
			ans_long_chain child_ans = getLongestChain(child, length+1);
			if((child_ans.length > max_length) || ((child_ans.length == max_length) && (child_ans.time < min_time)))
			{
				max_length = child_ans.length;
				min_time = child_ans.time;
				lastNode = child_ans.last_node;
			}
		}
		ans_long_chain final_ans(lastNode, min_time ,max_length);
		return final_ans;
	}
public:

	blockchain()
	{
		root = new node();
		root->blk->blockID = 0;
		root->parent = NULL;	
	}

	void add_block(int prevblockID, int blockID, double time_arrival)
	{
		node *lastNode = find_block(root, prevblockID);
		if(lastNode == NULL)
		{
			cout<<"Error"<<endl;
			return;
		}
		node *newNode = new node();
		newNode->parent = newNode;
		newNode->blk = new block();
		newNode->blk->prevblockID = prevblockID;
		newNode->blk->blockID = blockID;
		newNode->blk->time_arrival = time_arrival;
		lastNode->nextBlocks.push_back(newNode);
		// TODO : Manage list of transactions coming, adjust time based on latency..
	}

	void generate_node(double time)
	{
		node *current = new node();
		current->blk = new block();
		current->blk->blockID = block_counter;
		node *longestChain = (getLongestChain(root, 0).last_node);
		current->blk->prevblockID = longestChain->blk->blockID;
		current->blk->time_arrival = time;
		longestChain->nextBlocks.push_back(current);
		// we get the previous blockid
		block_counter++;
		// Now we need to broadcast this also 
	}
};	

void sorted_add(tnx trans, vector <tnx> &globalQueueTnx)
{
	if(globalQueueTnx.size() == 0)
		globalQueueTnx.push_back(trans);
	int i = 0;
	for( ; i < globalQueueTnx.size(); i++)
	{
		if(globalQueueTnx[i].time > trans.time)
			break;
	}
	globalQueueTnx.insert(globalQueueTnx.begin()+i, trans);
}

class network;
// Classes
class peer {
  private:
	int ID;
	// float amount;
	bool type;
	int activation;
	int connected;
	network *coin;
	vector <tnx> globalQueueTnx;
	vector<int> blocks_rec;
	bool anyBlockArrival = false; // set to false while generating time for next block while set true on receiving any block
	blockchain *chain;
  public:

	peer(int id, bool speed, int active, network *tmp) {
		ID = id;
		amount = 0;
		type = speed;
		activation = active;
		connected = 0;
		this->coin = tmp;
		chain = new blockchain();
	}

	void generate_transaction(double time) // for intial time would be zero
	{
		float pay_amt = amount * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		int rec = rand() % N;
		while(rec != ID)
		{
			rec = rand() % N;
		}
		tnx trans(++txn_counter, ID, rec, pay_amt, time);
		double next_time = time + exp_dist(txn_mean);
		event nextTrans(0, this->ID, next_time);
		time_simulator.sorted_event_add(nextTrans);
		// sorted_add(trans, globalQueueTnx);
		broadcast_tnx(trans, this, -1);
		// globalQueueTnx.add(txn_counter);
		// TODO : I have to broadcast it as well (Important)
	}
	void generate_block(double time)
	{

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

	void add_txn(tnx trans) {
		// need to change the time
		if(trans.send_id != this->ID)
		{
			double latency = coin->get_latency(mycoin->nodelist[trans.send_id], this, 0);
			trans.time += latency;

		}	
		sorted_add(trans, globalQueueTnx);
	}
	
	
	void add_blk(block blk)
	{
		if(blk.peer_id != this->ID)
		{
			double latency = coin->get_latency(mycoin->nodelist[blk.peer_id.send_id], this, 1); //TODO :  what should be sent instead of 1
			blk.time_arrival += latency;
		}
		blockchain.add_block(blk.prevblockID, blk.blockID, blk.time_arrival);
		// Handle Trans
	}

	bool txn_exists(int id) {
		for(auto it : globalQueueTnx)
			if(it.id == id)
				return true;
		return false;	
	}
	bool blk_exist(int id)
	{
		for(auto it : blocks_rec)
			if(it == id)
				return true;
		return false;	
	}
};


// Maintains the cryptocurrency network
class network {
  private:
  	
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

  	vector<peer*> nodelist;
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
			peer* newnode = new peer(id, type, active, this);
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

	// bool validate(tnx transaction) {
	// 	int send_id = transaction.send_id;
	// 	int txn_amount = transaction.amount;
	// 	peer* send_node = findnode(send_id);
	// 	if (send_node->get_amount() < txn_amount) 
	// 		return false;
	// 	else
	// 		return true;
	// }

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

	void broadcast_tnx(tnx &trans, peer* recv_node, int send_id) {
		if (recv_node->txn_exists(trans.id))
			return;
		else {
			recv_node->add_txn(trans);
			for( const auto& peerlist : adjlist[recv_node->get_id()] ) {
				if(peerlist.first->get_id() != send_id)
					broadcast_tnx(trans, peerlist.first, recv_node->get_id());
			}
		}
	}
	void broadcast_blk(block &blk, peer* recv_node, int send_id) {
		if (recv_node->blk_exists(blk.blockID))
			return;
		else {
			recv_node->add_blk(blk);
			for( const auto& peerlist : adjlist[recv_node->get_id()] ) {
				if(peerlist.first->get_id() != send_id)
					broadcast_blk(blk, peerlist.first, recv_node->get_id());
			}
		}
	}	
	// TODO : need to make a new function to broadcast blocks
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

	// Initialising the queue
	for(int i = 0; i < N; i++)
	{
		event txn(0, i, exp_dist(txn_mean));
		time_simulator.sorted_event_add(txn);
		event blk(1, i, exp_dist(block_mean));
		time_simulator.sorted_event_add(blk);
	}

	while(true)
	{
		if(time_simulator.size() != 0)
		{	
			event current = time_simulator[0];
			time_simulator.erase(time_simulator.begin());
			if(current.type == 0)
				mycoin->nodelist[current.peer_id].generate_transaction(current.time);
			else 
				mycoin->nodelist[current.peer_id].generate_block(current.time);	
  		}	  
	}
	// mycoin.print_graph();

	// cout << mycoin.get_latency(mycoin.findnode(0), mycoin.findnode(2), 1);

	// mycoin.broadcast(10, mycoin.findnode(0), 1);

	cout<<"\nend"<<endl;
	return 0;
}

// TODO : Should we make a block if low number of unspect trans
