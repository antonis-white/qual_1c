#include <vector>
#include <stdint.h>
#include <iostream>
#include <memory>

namespace subway {

enum sssp_algo_t { A_star_algo, spfa_algo };

using time_type = uint32_t; 
static inline constexpr size_t station_npos = -1;

struct line_station_t {
	size_t line; 
	size_t index;
	               
	line_station_t() = default;
	line_station_t(size_t line, size_t index) : line(line), index(index) {}
};
	
struct position_t : public line_station_t {
	size_t station_id;
	bool rev_order;
		
	position_t() = default;
	position_t(size_t line, size_t index, size_t station_id, bool rev_order)
		: line_station_t(line, index), station_id(station_id), rev_order(rev_order) {}
};

//pair-like data type
struct line_t {
	std::vector<size_t> order;
	std::vector<time_type> dist;
	time_type freq;

private:
	std::vector<time_type> cumsum;
public:
	line_t(const std::vector<size_t>& order, const std::vector<time_type>& dist, time_type freq) 
			: order(order), dist(dist), freq(freq), cumsum(order.size()) {

		if (dist.size() + 1 != order.size())
			throw std::logic_error("Wrong size of list of distances between stations in line");
		if (freq == 0)
			throw std::logic_error("Frequency of trains should be non-zero");	
		for (size_t i = 0; i < dist.size(); ++i) {
			if (dist[i] == 0)
				throw std::logic_error("Length of path should be non-zero");
			cumsum[i + 1] = cumsum[i] + dist[i];
		}
	}

	time_type closest_train(size_t index, bool rev_order, time_type cur_time) {
		time_type first_train = cumsum[index];
		if (rev_order)
			first_train = cumsum.back() - cumsum[index];
		if (first_train >= cur_time)
			return first_train;
		time_type rem = (cur_time - first_train) % freq;
		return cur_time + (rem == 0 ? 0 : freq - rem);
	}
};

/*
struct path_type {
	enum event_t { transfer, movement, start};

	struct event_type {
		event_t event;
		time_type end_time;
		position_t pos;

		event_type(event_t event, time_type end_time, position_t pos) :
			event(event), end_time(end_time), pos(pos) {}

		friend std::ostream& operator<<(std::ostream& out, const event_type& ev) {
			if (ev.event == event_t::transfer) {
				out << "Transfer to:\n";
			} else if (ev.event == event_t::movement) {
				out << "Movement to:\n";
			} else {
				out << "Starts from:\n";
			}
			out << "\tLine: " << ev.pos.line << '\n';
			out << "\tIndex: " << ev.pos.index << '\n';
			out << "\tStation: " << ev.pos.station_id << '\n';
			out << "\tGoing " << (ev.pos.rev_order ? "forward" : "backwards") << '\n';
			out << "\tCurrent time: " << ev.end_time;
			return out;
		}
	};
public:
	std::vector<event_type> events;
public:

	friend std::ostream& operator<<(std::ostream& out, const path_type& pth) {
		for (const auto& event : pth.events)
			out << event << '\n';
	    out << "Finished\n";
		return out;
	}
};
*/

using path_type = time_type;
	

template <sssp_algo_t sssp_algo = spfa_algo>
class subway_predictor {
private:
	std::vector<std::vector<line_station_t>> transfers_;
	std::vector<line_t> lines_;

public:
	size_t station_count() const { return transfers_.size(); }
	size_t line_count() const { return lines_.size(); }

private:
	struct sssp_state : public position_t {
		size_t trans_count;
		sssp_state() = default;
		sssp_state(size_t line, size_t index, size_t station_id, bool rev_order, size_t trans_cnt):
			position_t(line, index, station_id, rev_order), trans_count(trans_cnt) {} 
	};

	struct spfa_data_T {
		std::vector<sssp_state> q;
		std::vector<std::vector<std::vector<std::vector<time_type>>>> dist;
		std::vector<std::vector<std::vector<std::vector<sssp_state>>>> par;
		std::vector<std::vector<std::vector<std::vector<char>>>> in_q;	
	};
	struct A_star_data_T {};

	std::conditional_t<sssp_algo == spfa_algo, spfa_data_T, void*> spfa_data;
    std::conditional_t<sssp_algo == A_star_algo, A_star_data_T, void*> A_star_data;	

	//[line_id][index][ord][trans_cnt (lazy allocation)]
	template <typename cont_t> 
	void allocate_state_info(cont_t& cont) {
		cont.resize(line_count());
		for (size_t line_id = 0; line_id < line_count(); ++line_id) {
			cont[line_id].resize(lines_[line_id].order.size());
			for (auto& tmp : cont[line_id]) {
				tmp.resize(2); 
			}	
		}
	}

	size_t get_station_id(size_t line_id, size_t index) {
		return lines_[line_id].order[index];
	}

    template <bool consider_transfers>
	path_type A_star_impl(size_t from, size_t to, time_type start_time, time_type block) {
		path_type res;
		//todo
		throw std::runtime_error("Not implemented");
		return res; 
	}
 	
	template <bool consider_transfers>
	path_type spfa_impl(size_t from, size_t to, time_type start_time, time_type block) {   
		size_t delta_trans = consider_transfers;
		auto& q = spfa_data.q;
		auto& dist = spfa_data.dist;
		auto& par = spfa_data.par;
		auto& in_q = spfa_data.in_q;

		static bool was_allocated = false;
		if (!was_allocated) {
			allocate_state_info(dist);
			allocate_state_info(in_q);
			allocate_state_info(par);
			was_allocated = true;
		}

		size_t lit = 0;
		q.clear();

		auto relax = [&dist, &par, &q, &in_q, &block, this](
					size_t line, size_t ind, bool rev, size_t trans_cnt, time_type cur_time, const sssp_state& pstate) {
			if (cur_time > block)
				return;
			if (par[line][ind][rev].size() <= trans_cnt) {
		    	dist[line][ind][rev].resize(trans_cnt + 1, -1);
		    	in_q[line][ind][rev].resize(trans_cnt + 1, 0);
		    	par[line][ind][rev].resize(trans_cnt + 1);
		   	}
			if (dist[line][ind][rev][cur_time] > cur_time) {
				dist[line][ind][rev][trans_cnt] = cur_time;
				par[line][ind][rev][trans_cnt] = pstate;
				if (!in_q[line][ind][rev][trans_cnt])
					in_q[line][ind][rev][trans_cnt] = true;
				q.emplace_back(line, ind, this->get_station_id(line, ind), rev, trans_cnt);
			}
		};
		sssp_state dummy(-1, -1, -1, -1, -1);

		for (const auto& s: transfers_[from]) {
			for (int rev = 0; rev < 2; ++rev)
				relax(s.line, s.index, rev, 0, start_time, dummy);	
		}
		while (lit < q.size()) {
			auto& st = q[lit++];
			size_t station_id = get_station_id(st.line, st.index);
			time_type cur_time = dist[st.line][st.index][st.rev_order][st.trans_count];
			//transfers
			for (auto s: transfers_[station_id]) {
				for (int rev = 0; rev < 2; ++rev) {
					relax(s.line, s.index, rev, st.trans_count + delta_trans, cur_time, st); 
				}
			}
			auto nst = st;
			if (!st.rev_order && st.index + 1 == lines_[st.line].order.size())
				continue;
			if (st.rev_order && st.index == 0)
				continue;
			nst.index += (st.rev_order ? -1 : 1);
			cur_time = lines_[st.line].closest_train(st.index, st.rev_order, cur_time);
			cur_time += lines_[st.line].dist[st.rev_order ? nst.index : st.index];						
			relax(nst.line, nst.index, nst.rev_order, nst.trans_count, cur_time, st);
		}

		time_type ans = -1;
		for (const auto s: transfers_[to]) {
			for (int rev = 0; rev < 1; ++rev) {
				for (const auto& tmp : dist[s.line][s.index][rev])
					if (tmp < ans)
						ans = tmp;
			}
		}
		return ans - start_time;
	}

	template <bool consider_transfers>
	path_type sssp_impl(size_t from, size_t to, time_type start_time, time_type block) {
		if constexpr (sssp_algo == spfa_algo) {
			path_type ans = spfa_impl<consider_transfers>(from, to, start_time, block);
			return ans;
		}
		if constexpr (sssp_algo == A_star_algo) {
			path_type ans = A_star_impl<consider_transfers>(from, to, start_time, block);
			return ans;
		}
	}

public:
	subway_predictor(size_t station_cnt, const std::vector<line_t>& lines)
			: transfers_(station_cnt), lines_(lines) {
		for (size_t line_id = 0; line_id < line_count(); ++line_id) {
			const auto& order = lines_[line_id].order;
			for (size_t i = 0; i < order.size(); ++i) {
				size_t st = order[i];
				if (st >= station_count()) {
					throw std::logic_error("Wrong station id");
				}
				transfers_[st].emplace_back(line_id, i);
			}
		}
	}
	template <bool want_to_read=false>
	path_type predict(size_t from, size_t to, time_type start_time) { 
		path_type res = sssp_impl<false>(from, to, start_time, -1);
		if constexpr (!want_to_read)
			return res;
		//time_type spent_time = res.events.back().end_time - start_time; 
		time_type spent_time = res; 
		path_type res2 = sssp_impl<true>(from, to, start_time, start_time + 3 * spent_time / 2);
		return res2;
	}
}; 

} //namespace subway
