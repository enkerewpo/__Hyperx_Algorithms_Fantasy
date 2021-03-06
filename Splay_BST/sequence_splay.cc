/**
 *   sequence_splay.cc
 *
 * Copyright(C) 2017 Kvar_ispw17 project __Hyperx All rights reserved
 */

#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cassert>

const int N = 100 + 5;
int n, a[N];

class node;
node *root = NULL;
void tree_print(node*, int);
void print(node*);
node* search(node*, int);

class node {
public:
	/**
	 *   v    	the key value of the node  (int)
	 *   s     	the subtree size  (int)
	 *   rev    the tag of action 'reverse'  (bool)
	 *   c[]	the left and right child's pointers  (node*)
	 *   f 		the father's pointer  (node*)
	 */
	int v, s;
	bool rev;
	node *c[2], *f;
	node() { rev = 0, v = 0, s = 1, f = c[0] = c[1] = NULL; }
	node(int v, node *f) : v(v), f(f) { c[0] = c[1] = NULL, rev = 0, s = 1; }
	void up() {
		s = 1;
		if (c[0]) s += c[0]->s;
		if (c[1]) s += c[1]->s;
	}
	void down() {
		if (!rev) return;
		if (c[0]) c[0]->rev ^= 1;
		if (c[1]) c[1]->rev ^= 1;
		rev = 0;
		std::swap(c[0], c[1]);
	}
	void print() {
		printf("%d(%d) %d_%d |%d\n", v, s, c[0] ? c[0]->v : -1, c[1] ? c[1]->v : -1, f ? f->v : -1);
	}
};
/**
* access the relatoinship between son-node and father-node.
* @param  up   father-node
* @param  down son-node
* @return 	[if] son-node is the _right_ son
*/
int son(node* up, node* down) {	return down == up->c[1]; }
/**
* the element which indexed at position $k
* @param  k the $k position
* @return   the pointer
*/
node* at(int k) { return search(root, k); }
/**
* build sequence-splay from an array
* @param  o the root address of the generated splay
* @param  x the array
* @param  l left-bound closed
* @param  r right-bound closed
* @return   the pointer of the current sub-root
*/
node* build(node* &o, int x[], int l, int r) {
	if (l > r) return NULL;
	o = new node();
	if (l == r) { o->v = x[l]; return o; }
	int mid = ceil((l + r) / 2.);
	o->v = x[mid];
	o->c[0] = build(o->c[0], x, l, mid - 1); if (o->c[0]) o->c[0]->f = o;
	o->c[1] = build(o->c[1], x, mid + 1, r); if (o->c[1]) o->c[1]->f = o;
	o->up();
	return o;
}
/**
* this just will do what you think right now
* @param o pivot
* rotate pivot to its father's position.
*/
void rotate(node* o) {
	node* f = o->f, *gf = f->f;
	if (f == NULL) return;
	int d = son(f, o), gd = 0;
	if (gf) gd = son(gf, f);
	f->c[d] = o->c[d ^ 1]; if (o->c[d ^ 1]) o->c[d ^ 1]->f = f;
	o->c[d ^ 1] = f, f->f = o, gf ? gf->c[gd] = o : root = o, o->f = gf;
	o->up(), f->up();
}
/*	And you should never use the below one if you want to reverse your sequence.  */
void bad_rotate(node* &o, int d) {
	printf("rotate(%d, %s)\n", o->v, d ? "zag" : "zig");
	if (o->rev) {
		o->down();
		d ^= 1;
		puts("!");
	}
	node *f = o->f, *k = o->c[!d];
	k->down();
	o->c[!d] = k->c[d];
	if (o->c[!d]) o->c[!d]->f = o;
	k->c[d] = o, k->f = f, o->f = k, f ? f->c[son(f, o)] = k : root = k; o->up(), k->up(), o = k;
}
/**
* search the kth element and push down the tag on the access path
* @param  o current visiting node
* @param  k current goal
* @return   the pointer
*/
node* search(node* o, int k) {
	o->down();
	int s = 0;
	if (o->c[0]) s = o->c[0]->s;
	if (k < s + 1) return search(o->c[0], k);
	else if (k > s + 1) return search(o->c[1], k - (s + 1));
	else return o;
}
/**
* splay anything to anywhere
* @param o node needed splayed
* @param t node you want $o to splay 'under'
* after the splay, $o will be a son of $t.
*/
void splay(node* o, node* t) {
	while (o->f != t) {
		node* f = o->f;
		node* gf = f->f;
		if (gf == t) { rotate(o); break; }
		int d = son(f, o);
		int gd = son(gf, f);
		d == gd ? rotate(f) : rotate(o);
		rotate(o);
	}
}
/**
* insert a externel sequence 'into' the splay
* @param o host
* @param k join exactly after the kth element
* @param x array
* @param l left-bound closed
* @param r right-bound closed
*/
void insert(node* o, int k, int x[], int l, int r) {
	node* t = build(t, x, l, r), *p = at(k);
	splay(p, NULL);	node* q = p->c[1]; p->c[1] = t, t->f = p;
	splay(p = at(k + t->s), NULL); p->c[1] = q, q->f = p;
}
/**
* erase the sequence $[l, r] on splay
* @param o host
* @param l left-bound closed
* @param r right-bound closed
*/
void erase(node* o, int l, int r) {
	node* p = at(l - 1), *q = at(r + 1);
	splay(p, NULL);
	splay(q, p);
	q->c[0] = NULL;
}
/**
 * reverse the sequence $[l, r]
 * @param o host
 * @param l left-bound closed
 * @param r right-bound closed
 */
void reverse(node* o, int l, int r) {
	node *p = at(l - 1), *q = at(r + 1);
	splay(p, NULL);
	splay(q, p);
	q->c[0]->rev ^= 1;
}
/**
 * print by 'middle-first' order
 * @param o node
 */
void print(node* o) {
	o->down();
	if (o->c[0]) print(o->c[0]);
	printf("%d ", o->v);
	if (o->c[1]) print(o->c[1]);
}
/**
 * print by 'left-first' order
 * @param o   node
 * @param dep depth
 */
void tree_print(node* o, int dep) {
	for (int i = 0; i < dep; i++) printf("    ");
	o->print();
	if (o->c[0]) tree_print(o->c[0], dep + 1);
	if (o->c[1]) tree_print(o->c[1], dep + 1);
}

/**
 * [ follow the guide in main() to get to know this program ]
 */
int main() {
	/*   Build the splay from the $n elements in array $a[]   */
	printf("Input the numbers count of your sequence: "), scanf("%d", &n);
	for (int i = 1; i <= n; i++) printf("[%d] = ", i), scanf("%d", &a[i]);
	build(root, a, 0, n + 1);
	printf("The extended sequence: "), print(root), printf("\n");

	/*   Insert $m elements in array $b[] after the $k th element in the splay hosted by $root   */
	int b[N], k, m;
	printf("Input the numbers count of the sequence for inserting: "), scanf("%d", &m);
	for (int i = 1; i <= m; i++) printf("[%d] = ", i), scanf("%d", &b[i]);
	printf("The position you want to be the start of the inserting: (ignore the '0') "), scanf("%d", &k);
	insert(root, k, b, 1, m);
	printf("The sequence after insert: "), print(root), printf("\n");

	/*  Erase the sequence $[eL, eR] in $root  */
	int eL, eR;
	printf("Input the left bound and the right bound for your erasing: (ignore the '0')"), scanf("%d%d", &eL, &eR);
	erase(root, eL + 1, eR + 1);
	printf("The sequence after erase: "), print(root), printf("\n");

	/*  Reverse the sequence T time between [L_i, R_i]  */
	int T;
	printf("Input how many times you want to reverse: "), scanf("%d", &T);
	while (T--) {
		int L, R;
		printf("[%d times left] ", T + 1);
		printf("Input the left and right bound for current reverse action: (ignore the '0')"), scanf("%d%d", &L, &R);
		reverse(root, L + 1, R + 1);
	}
	printf("The final sequence: "), print(root), printf("\n");
	return 0;
}
