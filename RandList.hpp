// #pragma once
#include <algorithm>
#include <functional>
#include <vector>
#include <unordered_map>
//using std::udi;
using std::vector;
using udi = unsigned int;
class RandList {
   public:
    udi cap;
    udi vnum;
    vector<udi> vlist;
    vector<udi> vpos;
    //unordered_map<udi,udi> syn;
    RandList(udi _cap = 0) : cap(_cap), vnum(0) {
        vlist.reserve(cap);
        //syn.reserve(cap);
        vpos.resize(cap);
        std::fill(vpos.begin(), vpos.end(), cap);
    }
    void swap(RandList& other) noexcept {
        using std::swap;
        swap(cap, other.cap);
        swap(vnum, other.vnum);
        swap(vlist, other.vlist);
        swap(vpos, other.vpos);
        //syn.swap(other.syn);
    }
    inline void init(udi _cap){
        cap=_cap;
        vnum=0;
        vlist.reserve(cap);
        vpos.resize(cap);
        std::fill(vpos.begin(), vpos.end(), cap);
        //syn.reserve(cap);
    }
    RandList& operator=(const RandList& other) {
        if (this != &other) {
            cap = other.cap;
            vnum = other.vnum;
            vlist = other.vlist;
            vpos=other.vpos;
        }
        return *this;
    }
    inline bool contains(udi vid) const { 
        if(vid>cap) return false;
        return vpos[vid] != cap; 
        }
    inline void add(udi vid) {
        if(vid>cap) return;
        if(contains(vid)) return;
        vlist.emplace_back(vid);
        vpos[vid] = vnum++;
    }
    inline void remove(udi vid) {
        if(vid>cap) return;
        if(!contains(vid)) return;
        udi last_id = vlist[vnum - 1];
        udi id_pos = vpos[vid];
        vlist[id_pos] = last_id;
        vpos[last_id] = id_pos;

        vlist.pop_back();
        --vnum;
        vpos[vid] = cap;
    }
    inline void clear() {
        for (udi i = 0; i < vnum; i++) {
            vpos[vlist[i]] = cap;
        }
        vlist.clear();
        vnum = 0;
    }

    void valueEqual(const RandList& other) {
        clear();
        for (udi i = 0; i < other.vnum; i++) {
            udi vid = other.vlist[i];
            add(vid);
        }
    }

    void reinit(udi _cap) {
        cap = _cap;
        clear();
        vlist.reserve(cap);
        vpos.resize(cap);
        std::fill(vpos.begin(), vpos.end(), cap);
    }
    inline bool empty() const {
        return vnum == 0;
    }
};
