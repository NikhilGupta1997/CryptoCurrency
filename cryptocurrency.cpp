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
using namespace std;

// Global Variables
int N = 10; //Number of peers
float txn_mean = 5.0; // Mean of exponential transaction distribution function
float block_mean = 10.0; // TODO : Nikhil
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
	vector<tnx> spent;
};
struct event
{
	int event_type = 0; //0 for tranx, 1 for block generation, 2 for add block 
	int peer_id;
	double time;
	double time_gen = 0.0; 
	block * blk; 
	event(int event_type, int peer_id, double time)
	{
		this->event_type = event_type;
		this->peer_id = peer_id;
		this->time = time;
		blk = NULL;
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

struct node
{
	block *blk;
	node *parent;
	vector<node *> nextBlocks;
	vector<float> peer_amount;
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
	public:
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

	node* find_common_anc(node *prev, node *curr)
	{
		node *tmp = prev;
		unordered_set<node*> parents;
		while(tmp !=NULL)
		{
			// cout<<"yo : "<< tmp->blk->blockID<<endl;
			parents.insert(tmp->parent);
			tmp = tmp->parent;
		}
		node* tempo = curr;
		while(tempo)
		{
			if(parents.find(tempo) != parents.end())
				return tempo;
			tempo = tempo->parent; 
		}
		return NULL; // Never gonna happen
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

	blockchain()
	{
		root = new node();
		root->blk = new block();
		root->blk->blockID = 0;
		root->parent = NULL;	
	}

	void update(vector <tnx> &globalQueueTnx, node* prev, node* curr, int peer_id)
	{
		// cout<<"update"<<endl;
		node *anc = find_common_anc(prev, curr);
		// cout<<"got anc";
		node *tmp = prev;
		while(tmp != anc)
		{
			// i am adding all unspent here which is wrong but ta ka kat denge
			for(auto it : tmp->blk->spent)
				sorted_add(it, globalQueueTnx);
			tmp = tmp->parent;
		}
		tmp = curr;
		while(tmp != anc)
		{
			unordered_set<int> trans; 
			for(auto it : tmp->blk->spent)
				trans.insert(it.id);
			for(int i = 0; i < globalQueueTnx.size(); i++)
				if(trans.find(globalQueueTnx[i].id) != trans.end())
					globalQueueTnx.erase(globalQueueTnx.begin()+i);
			tmp = tmp->parent; 
		}
		// cout<<"update done"<<endl;
	}

	node *add_block(int prevblockID, int blockID, double time_arrival)
	{
		node *lastNode = find_block(root, prevblockID);
		if(lastNode == NULL)
		{
			cout<<"Error : Block id "<<prevblockID<<endl;
			return NULL;
		}
		node *newNode = new node();
		newNode->parent = lastNode;
		newNode->blk = new block();
		newNode->blk->prevblockID = prevblockID;
		newNode->blk->blockID = blockID;
		newNode->blk->time_arrival = time_arrival;
		lastNode->nextBlocks.push_back(newNode);
		return newNode;
	}

	node * generate_node(double time)
	{
		node *current = new node();
		current->blk = new block();
		current->blk->blockID = block_counter;
		ans_long_chain long_chain = getLongestChain(root, 0);
		node *longestChain = (long_chain.last_node);
		current->parent = longestChain;
		current->blk->prevblockID = longestChain->blk->blockID;
		// cout<<"Block no is "<<block_counter<<" prev was "<< current->blk->prevblockID<<endl;
		current->blk->time_arrival = time;
		longestChain->nextBlocks.push_back(current);
		// we get the previous blockid
		block_counter++;
		return current;
	}
};	


// TODO : Remove Forward dependency : 
// https://stackoverflow.com/questions/19962812/error-member-access-into-incomplete-type-forward-declaration-of
class network;
// Classes
class peer {
  private:
	int ID;
	// float amount;
	bool type;
	int activation;
	int connected;
	// network *coin;
	vector <tnx> globalQueueTnx; // to be maintained at the current block
	vector<int> blocks_rec;
	// unordered_set<int> orphans;
	unordered_map<int, vector<block>> orphan_blk;
	double lastBlockArrival = 0.0; 
	blockchain *chain;
	node *longest_one;
  public:

	peer(int id, bool speed, int active) {
		ID = id;
		// amount = 0;
		type = speed;
		activation = active;
		connected = 0;
		chain = new blockchain();
		longest_one = chain->root;
	}
	// network* give_network();
	// void set_network(network *coin);
	
	void generate_transaction(double time, network * tmp); // for intial time would be zero

	void generate_block(double time, double time_gen, network * tmp);

	int get_id() { // Get node id
		return ID;
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

	void add_txn(tnx trans, network * tmp); 
		
	void add_blk(block &blk);

	void add_blk_sim(block &blk, network * tmp);

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
			peer* newnode = new peer(id, type, active);
			// newnode->set_network(this);
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
		// cout << "\np = " << p << endl;
		// cout << "c = " << c << endl;
		// cout << "d_mean = " << d_mean << endl;
		float total_latency = (float)m_size / c + exp_dist(d_mean) + p;
		return total_latency;
	}

	void broadcast_tnx(tnx &trans, peer* recv_node, int send_id) {
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
	void broadcast_blk(block &blk, peer* recv_node, int send_id) {
		if (recv_node->blk_exist(blk.blockID))
			return;
		else {
			recv_node->add_blk_sim(blk, this);
			for( const auto& peerlist : adjlist[recv_node->get_id()] ) {
				if(peerlist.first->get_id() != send_id)
					broadcast_blk(blk, peerlist.first, recv_node->get_id());
			}
		}
	}	
	
};

void peer::generate_transaction(double time, network * tmp) // for intial time would be zero
{

	float pay_amt = longest_one->peer_amount[this->ID] * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
	int rec = rand() % N;
	while(rec == ID)
	{
		rec = rand() % N;
	}
	tnx trans(++txn_counter, ID, rec, pay_amt, time);
	double next_time = time + exp_dist(txn_mean);
	event nextTrans(0, this->ID, next_time);
	sorted_event_add(nextTrans);
	tmp->broadcast_tnx(trans, this, -1);
}
void peer::generate_block(double time, double time_gen, network * tmp)
{
	// add only those transaction for which time is less than this
	// cout<<"Generating Block"<<endl;
	if(lastBlockArrival <= time && lastBlockArrival > time_gen)
		return;
	if(globalQueueTnx.size() == 0)
		return;
	// cout<<"Pass"<<endl;
	node *current = chain->generate_node(time);
	current->blk->peer_id = this->ID;
	int i;
	for(i = 0; i < globalQueueTnx.size(); i++)
	{
		if(globalQueueTnx[i].time <= time)
			current->blk->spent.push_back(globalQueueTnx[i]);
		else 
			break;
	}
	globalQueueTnx.erase(globalQueueTnx.begin(), globalQueueTnx.begin()+i);
	// cout<<"Block No "<<block_counter<<" Gen by "<<this->ID<<" with transactions = "<<current->blk->spent.size()<<endl;
	tmp->broadcast_blk(*(current->blk), this, -1);
}
void peer::add_txn(tnx trans, network *tmp) {
	// need to change the time
	if(trans.send_id != this->ID)
	{
		double latency = tmp->get_latency(tmp->nodelist[trans.send_id], this, 0);
		trans.time += latency;
	}	
	sorted_add(trans, globalQueueTnx);
}
void peer::add_blk_sim(block &blk, network *tmp)
{
	blocks_rec.push_back(blk.blockID);
	if(blk.peer_id != this->ID)
	{
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
void peer::add_blk(block &blk)
{
	node *current = chain->add_block(blk.prevblockID, blk.blockID, blk.time_arrival);
	bool flag = false;
	if(!current) 
	{
		// orphans.insert(blk.prevblockID);
		orphan_blk[blk.prevblockID].push_back(blk);
		flag = true;
	}	
	if(blk.peer_id != this->ID)
	{
		// cout<<"1"<<endl;
		lastBlockArrival = blk.time_arrival;
		// generate a random variable and create event for next block gen
		event blk_gen(1, this->ID, lastBlockArrival+exp_dist(block_mean));
		blk_gen.time_gen = lastBlockArrival;
		sorted_event_add(blk_gen);
		if(flag)
			return;
		// cout<<"2"<<endl;
		unordered_set<int> trans_id;
		for(auto it : blk.spent)
			trans_id.insert(it.id);
		// cout<<"3"<<endl;
		for(int i = 0; i < globalQueueTnx.size(); i++)
		{
			if(trans_id.find(globalQueueTnx[i].id) != trans_id.end())
			{
				globalQueueTnx.erase(globalQueueTnx.begin() + i);
				i--;
			}
		}
	}
	if(flag)
		return;
	// cout<<"4"<<endl;

	// cout<<"5"<<endl;
	ans_long_chain long_chain = chain->getLongestChain(chain->root, 0);
	// cout<<"6"<<endl;
	node* new_longest_one = long_chain.last_node;
	if(this->ID == 0)
	{	cout<<"For peer "<<this->ID<< ", length is "<<long_chain.length<<endl;
		if(blk.prevblockID != longest_one->blk->blockID)
			cout<<"Forking"<<endl;
	}
	if(blk.prevblockID != longest_one->blk->blockID && new_longest_one->blk->blockID == blk.blockID)	
		chain->update(globalQueueTnx, longest_one, new_longest_one, this->ID);
	longest_one = new_longest_one;
	// cout<<"7"<<endl;
	for(int i = 0; i < N; i++)
		current->peer_amount[i] = current->parent->peer_amount[i];
	for(auto it : blk.spent)
	{
		current->blk->spent.push_back(it);
		if(it.amount > current->peer_amount[it.send_id])
			continue;
		current->peer_amount[it.send_id] -= it.amount;
		current->peer_amount[it.recv_id] += it.amount;
		current->peer_amount[current->blk->peer_id] += 50;
	}
	if(orphan_blk.find(current->blk->blockID) != orphan_blk.end())
	{
		for(auto it: orphan_blk[current->blk->blockID])
		{
			// cout<<"Retry : block "<<it.blockID<<endl;
			this->add_blk(it);
		}
		orphan_blk.erase(current->blk->blockID);
	}
}
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
		sorted_event_add(txn);
		event blk(1, i, exp_dist(block_mean));
		sorted_event_add(blk);
	}
	int counter = 0 ;
	while(true)
	{
		// if(counter %500 == 0)
		// {

		// }
		counter++;
		if(time_simulator.size() != 0)
		{	
			// cout<<"Num_blks : "<<block_counter<<" , Num_tnx: "<<txn_counter<<endl;
			event current = time_simulator[0];
			// cout<<current.event_type<<endl;
			time_simulator.erase(time_simulator.begin());
			if(current.event_type == 0)
				mycoin.nodelist[current.peer_id]->generate_transaction(current.time, &mycoin);
			else if(current.event_type ==1)
				mycoin.nodelist[current.peer_id]->generate_block(current.time, current.time_gen, &mycoin);	
  			else if(current.event_type ==2)
				mycoin.nodelist[current.peer_id]->add_blk(*(current.blk));	
  		
  		}	  
	}
	// mycoin.print_graph();

	// cout << mycoin.get_latency(mycoin.findnode(0), mycoin.findnode(2), 1);

	// mycoin.broadcast(10, mycoin.findnode(0), 1);

	cout<<"\nend"<<endl;
	return 0;
}


