from collections import deque
from cStringIO import StringIO
from copy import deepcopy
import re

IS_FLOAT = re.compile('^\s*(?:-?((\d*\.)?\d+)|((\d*\.\d+|\d+)(e|E)(\+|-)?\d+))\s*$') 
token_pat = re.compile('(?:\(|\)|,|\:|;|[^():;,]+)')

class ParseType:

    TREE='TREE'
    BRANCH='BRANCH'
    LEAF='LEAF'
    STRING='STRING'
    FLOAT='FLOAT'


class Leaf(object):
    
    def __init__(self, label, length=None):
        self.label = label
        self.length = length
        self.ancestor = None

    def setAncestor(self, anc):
        self.ancestor = anc

    def getAncestor(self):
        return self.ancestor

    def __str__(self):
        return '%s%s' % (self.label, self.length != None and ':%.5f' %self.length or '')


class Branch(object):

    def __init__(self, subtrees, label=None, length=None):
        self.subtrees = list()
        self.ancestor = None
        for t in subtrees:
            if type(t) != Branch and type(t) != Leaf:
                raise TypeError('%s is not of type Leaf or Branch but of %s' %(t, type(t)))
            t.setAncestor(self)
            self.subtrees.append(t)
        self.label = label or ''
        self.length = length

    def setAncestor(self, anc):
        self.ancestor = anc

    def getAncestor(self):
        return self.ancestor

    def __str__(self):
        return '(%s)%s%s' % (','.join(map(str, self.subtrees)),self.label != None
                and self.label or '', self.length != None and ':%.5f' %self.length or '')

    def getLeaves(self):
        res = []
        for t in self.subtrees:
            if type(t) == Leaf:
                res.append(t)
            elif type(t) == Branch:
                res.extend(t.getLeaves())
        return res

    def getNodes(self):
        res = []
        for t in self.subtrees:
            res.append(t)
            if type(t) == Branch:
                res.extend(t.getNodes())
        return res


class Tree(object):
    def __init__(self, subtree):
        self.subtree = subtree
    def __str__(self):
        return str(self.subtree) + ';' 
    def getLeaves(self):
        return self.subtree.getLeaves()
    def getNodes(self):
        return self.subtree.getNodes()


def parse_tree_iterator(data):
    if type(data) == str:
        data = StringIO(data)
    cursor = data.read(1)
    program = StringIO()
    while cursor:
        program.write(cursor)
        next_c = data.read(1)
        if cursor==';' or not next_c:
            t = __parse_tree__(program.getvalue())
            if t:
                yield t[0]
            program = StringIO()
        cursor = next_c

def parse_tree(data):
    res = []
    if type(data) == str:
        data = StringIO(data)
    cursor = data.read(1)
    program = StringIO()
    while cursor:
        program.write(cursor)
        next_c = data.read(1)
        if cursor==';' or not next_c:
            res.extend(__parse_tree__(program.getvalue()))
            program = StringIO()
        cursor = next_c
    if len(res) == 1:
        return res[0]
    return res


def __parse_tree__(program):
    res = []
    stack = [ParseType.TREE]
    for token in token_pat.findall(program):
        prev = stack[-1]
        if token == '(':
            if prev != ParseType.LEAF and prev != ParseType.TREE:
                raise SyntaxError('could not parse tree')
            # remove faulty prediction
            if prev == ParseType.LEAF:
                stack.pop()
            stack.append(ParseType.BRANCH)
            stack.append(ParseType.LEAF)
        elif token == ')':
            prev = stack.pop()
            subtrees = deque()
            while(prev != ParseType.BRANCH):
                subtrees.appendleft(prev)
                if not len(stack):
                    raise SyntaxError('too many closing brackets..')
                prev = stack.pop()
            stack.append(Branch(subtrees))
        elif token == ',':
            # predict leaf
            stack.append(ParseType.LEAF)
        elif token == ':':
            stack.append(ParseType.FLOAT)
        elif token == ';':
            if len(stack) > 1:
                prev = stack.pop()
                subtrees = deque()
                while(prev != ParseType.TREE):
                    if type(prev) != Leaf and type(prev) != Branch:
                        raise SyntaxError('unexpected end of tree')
                        continue
                    subtrees.appendleft(prev)
                    prev = stack.pop()
                if len(subtrees) > 1:
                    subtrees = [Branch(subtrees)]
                stack.append(Tree(subtrees[0]))
            res.append(stack.pop())
            if len(stack) != 0:
                raise SyntaxError('could not reduce tree...')
            stack.append(ParseType.TREE)
        elif not token.strip():
            # ignore whitespaces
            continue
        elif prev == ParseType.LEAF:
            # remove LEAF flag
            stack.pop()
            stack.append(Leaf(token.strip()))
        elif prev == ParseType.FLOAT:
            if IS_FLOAT.match(token):
                # remove FLOAT flag
                stack.pop()
                stack[-1].length = float(token)
            else:
                raise SyntaxError('%s was expected to be a number' %token)
        elif type(prev) == Branch:
            prev.label = token.strip()
        elif prev == ParseType.TREE:
            stack.append(Leaf(token.strip()))
        else:
            raise SyntaxError('unknown operator')
    if len(stack) > 1 or (stack and stack[0] != ParseType.TREE):
        raise SyntaxError('Error while parsing tree')
    return res


def constructDistTable(tree):
    res = {}
    leaves = [(l, [(0, l.label)]) for l in tree.getLeaves()]
    while len(leaves) > 2:
        hasMatch = False
        for i in range(len(leaves)-1):
            for j in range(i+1, len(leaves)):
                if leaves[i][0].getAncestor() == leaves[j][0].getAncestor():
                    joinItems(leaves, i, j, res)
                    hasMatch = True
                    break;
            if hasMatch:
                break
    joinItems(leaves, 0, 1, res)
    return res


def joinItems(leaves, i, j, res):
    l1 = leaves.pop(j)
    l2 = leaves.pop(i)
    anc = l1[0].getAncestor()
    for d1, k in l1[1]:
        for d2, l in l2[1]:
            key = [k, l]
            key.sort()
            res[tuple(key)] = d1 + d2 + l1[0].length + l2[0].length
    k = [(d+l1[0].length, n) for d, n in l1[1]]
    l = [(d+l2[0].length, n) for d, n in l2[1]]
    leaves.append((anc, k+l))
                            

def calculateSpan(tree):
    dists = constructDistTable(tree.subtree)
    return max(dists.values())


def reroot(new_root):
    br = deepcopy(new_root)
    if type(br) == Leaf:
        tmp = Branch(list(), label=br.label, length=br.length)
        tmp2 = br.getAncestor()
        tmp2.subtrees.remove(br)
        tmp2.subtrees.append(tmp)
        tmp.setAncestor(tmp2)
        br = tmp

    res  = Tree(br)

    cur = br
    new_anc = None
    new_len = None

    while cur:
        anc = cur.getAncestor()
        cur.setAncestor(new_anc)
        if anc != None:
            cur.subtrees.append(anc)
            anc.subtrees.remove(cur)
            tmp = cur.length
            cur.length = new_len
            new_len = tmp
        elif len(cur.subtrees) == 1:
            single = cur.subtrees[0]
            new_anc.subtrees.remove(cur)
            new_anc.subtrees.append(single)
            single.setAncestor(new_anc)
            if new_len:
                single.length = (single.length != None and single.length or 0) \
                    + (new_len != None and new_len or 0)
        tmp = cur
        cur = anc
        new_anc = tmp

    return res


def rescale_absolute(tree, max_length):
    span = calculateSpan(tree)
    s = max_length/float(span)
    return rescale(tree, s)


def rescale(tree, scale_factor):

    res = deepcopy(tree)

    stack = list()
    stack.append(res.subtree)
    
    while stack:
        x = stack.pop()
        if x.length:
            x.length *= scale_factor
        if type(x) == Branch:
            stack.extend(x.subtrees)

    return res

