/*
 *
 * @brief This is a splay BST implement for editing the sequence.
 * @date 2017.12.10
 * @author Kvar_ispw17
 *
 *  balanced_splay.cc  Copyright (C) 2017 Kvar_ispw17
 *
 */

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <algorithm>

class node {
public:

	class data {
	public:
		int idx;
		int val;
		data() { idx = val = 0; }
		data(int x) { idx = 0, val = x; }
		data(int x, int y) { idx = x, val = y; }
		bool operator > (const data& rhs) const {
			return idx > rhs.idx;
		}
	};
	
	data key;      /* key value */
	int siz;       /* size of subtree */
	node* fa;      /* pointer to father */
	node* ch[2];   /* pointers to children nodes */

	node() { fa = ch[0] = ch[1] = NULL, key = 0, siz = 1; }
	node(data _key = 0, node* _fa = NULL) : key(_key), fa(_fa) { ch[0] = ch[1] = NULL, siz = 1; }
	void maintain() {
		siz = 1;
		if (ch[0] != NULL) siz += ch[0]->siz;
		if (ch[1] != NULL) siz += ch[1]->siz;
	}
};

node* root = NULL;

int son(node* f, node* s) {
	return f->ch[1] == s;
}

node* search(node* o, int k) {
	int size = o->siz;
	if (k < size) return search(o->ch[0], k);
	else if (k > size) return search(o->ch[1], k - size);
	else return o;
}

void rotate(node* &o, int d) {
	node* k = o->ch[d ^ 1];
	node* f = o->fa;
	o->ch[d ^ 1] = k->ch[d];
	if (o->ch[d ^ 1]) o->ch[d ^ 1]->fa = o;
	f == NULL ? root = k : f->ch[son(f, o)] = k;
	k->ch[d] = o;
	k->fa = f;
	o->fa = k;
	o = k;
	o->maintain();
	k->maintain();
}

void splay(node* o, node* t) {
	while (o->fa != t) {
		node* f = o->fa;
		node* gf = f->fa;
		if (gf == t) rotate(f, !son(f, o));
		else {
			if (son(gf, f) ^ son(f, o))
				rotate(f, !son(f, o)), rotate(gf, !son(gf, f));
			else rotate(gf, !son(gf, f)), rotate(f, !son(f, o));
		}
	}
}

void insert(node::data val) {
	if (root == NULL) {
		root = new node(val, NULL);
		return;
	}
	node* p;
	int d;
	for (p = root, d = val > p->key; ; p = p->ch[d], d = val > p->key) {
		if (p->ch[d] == NULL) {
			p->ch[d] = new node(val, p);
			splay(p->ch[d], NULL);
			return;
		}
	}
}

void insert(int x, int idx) {
	node::data t(idx, x);
	insert(t);
}

void print(node* o) {
	if (o->ch[0] != NULL) print(o->ch[0]);
	printf("%d ", o->key.val);
	if (o->ch[1] != NULL) print(o->ch[1]);
}

void build(int a[], int n) {
	for (int i = 1; i <= n; i++) insert(a[i], i);
}

void erase(int idx) {
	node* o = search(root, idx);
	if(o != NULL) {
		splay(o, NULL);
		if(o->ch[0] == NULL) {
			root = o->ch[1];
			if(root) root->fa = NULL;
			root->maintain();
		} else {
			node* p = o->ch[0];
			while(p->ch[1]) p = p->ch[1];
			splay(p, o);
			p->fa = NULL;
			root = p;
			p->ch[1] = o->ch[1];
			if(p->ch[1]) p->ch[1]->fa = p;
			p->maintain();
		}
	}
}

void split(int k, node* &nroot) {
	node* o = search(root, k);
	splay(o, NULL);
	nroot = o->ch[1];
	o->ch[1] = NULL;
	o->maintain();
}

void merge(node* lhs, node* rhs) {
	node* o = search(root, lhs->siz);
	splay(o, NULL);
	o->ch[1] = rhs;
	o->maintain();
}

int main() {
	// Build the sequence
	int n;
	scanf("%d", &n);
	int a[100];
	for (int i = 1; i <= n; i++) scanf("%d", &a[i]);
	build(a, n);
	print(root), printf("\n");

	// Erase the kth element
	int p;
	scanf("%d", &p);
	erase(p);
	print(root), printf("\n");
	
	// Split from the kth element
	int k;
	node* newroot = NULL;
	scanf("%d", &k);
	split(k, newroot);
	print(newroot), printf("\n");
	
	// Merge them back
	merge(root, newroot);
	print(root), printf("\n");
	return 0;
}
