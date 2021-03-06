/**
 *  link-cut_tree.cc
 *  Copyright(C) 2017 Kvar_ispw17, project __Hyperx All rights reserved
 */

#include <cstring>
#include <algorithm>
#include <cstdio>

#define null 0x0
const int N = 1000 + 5;

class node {
public:
    /**
     *   v          the value of the node
     *   s          the size of the subtree
     *   fa         the address of the parent node
     *   path_fa    the edge on origin tree
     *   ch[]       the addresses of the children node
     */
    int v, s, fa, path_fa, ch[2];
    bool rev;
    node() { ch[0] = ch[1] = fa = v = s = rev = null; }
    node(int v) : v(v) { ch[0] = ch[1] = fa = rev = null, s = 1; }
    node(int v, int fa) : v(v), fa(fa) { ch[0] = ch[1] = rev = null, s = 1; }
} tr[N];
 
void up(int u) {
    tr[u].s = 1;
    if(tr[u].ch[0]) tr[u].s += tr[tr[u].ch[0]].s;
    if(tr[u].ch[1]) tr[u].s += tr[tr[u].ch[1]].s;
}
void push_down(int u) {
    if(tr[u].rev) {
        if(tr[u].ch[0]) tr[tr[u].ch[0]].rev ^= 1;
        if(tr[u].ch[1]) tr[tr[u].ch[1]].rev ^= 1;
        tr[u].rev = 0;
        std::swap(tr[u].ch[0], tr[u].ch[1]);
    }
}
    
int son(int u) { return tr[tr[u].fa].ch[1] == u; }
int isroot(int u) { return tr[tr[u].fa].ch[0] != u && tr[tr[u].fa].ch[1] != u; }
void down(int u) { if(tr[u].fa != null) down(tr[u].fa); push_down(u); }
    
void rotate(int u) {
    int f = tr[u].fa, pf = tr[f].fa, d = son(u);
    tr[f].ch[d] = tr[u].ch[d ^ 1];
    if (tr[f].ch[d]) tr[tr[f].ch[d]].fa = f;
    tr[u].ch[d ^ 1] = f;
    tr[f].fa = u;
    tr[u].fa = pf;
    if (!isroot(f)) tr[pf].ch[son(f)] = u;
}
void splay(int u) {
    down(u);
    for(int f = tr[u].fa; !isroot(u); rotate(u))
        if (!isroot(f)) rotate(son(u) == son(f) ? f : u);
}

/* Make the path from the 'root' of u to u all become [preferred edge]. */
void access(int u) {
    splay(u);
    tr[u].ch[1] = null;
    while(tr[u].path_fa) {
        int f = tr[u].path_fa;
        splay(f);
        tr[f].ch[1] = u;
        tr[u].path_fa = null;
        tr[u].fa = f;
        u = f;
    }
}
    
void makeroot(int u) { access(u); splay(u); tr[u].rev ^= 1; }
int getroot(int u) {
    access(u);
    splay(u);
    while(tr[u].ch[0]) u = tr[u].ch[0];
    return splay(u), u;
}
    
void link(int u, int v) { makeroot(u); tr[u]->path_fa = v; }
void cut(int u, int v) {
    makeroot(u);
    access(v);
    splay(v);
    tr[u].fa = null;
    tr[v].ch[0] = null;
}

int main() {
	return 0;
}
