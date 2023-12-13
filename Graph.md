# Graph Algorithm
## Graph Traversal
### Depth First Search

<i>Pseudocode: </i>

```
DFS(node):
    mark node as visited

    /*
        process the node
    */

    for each neighbor of node:
        if neighbor is not visited:
            DFS(neighbor)
```

<b>DFS on adjacency list</b>

```c++
void DFS(int node, &vector<bool> visited, vector<int> adjList[]){
    visited[node] = true;
    for (int neighbor : adjList[node]){
        if (!visited[neighbor])
            DFS(neighbor, visited, adjList);
    }
}
```

<b>DFS on adjacency matrix</b>

```c++
void DFS(int node, &vector<bool> visited, int adjMat[N][N]){
    visited[node] = true;
    for (int neighbor=0;neighbor<N;neighbor++){
        if (!visited[neighbor] && adjMat[node][neighbor] != 0)
            DFS(neighbor, visited, adjMat);
    }
}
```

<b>DFS on grid</b>

```c++
//                               Up     Down     Left    Right
vector<pair<int,int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

void DFS(int curr_row, int curr_col, bool visited[N][N], int grid[N][N]){
    visited[curr_row][curr_col] = true;
    for (auto [add_row, add_col] : directions){
        int next_row = curr_row + add_row;
        int next_col = curr_col + add_col;
        if (next_row < 0 || next_col < 0 || next_row >= ROW || next_col >= COL || visited[next_row][next_col])
            continue;
        DFS(next_row, next_col, visited, grid);
    }
}
```


### Breath First Search

<i>Pseudocode:</i>

```
initialize an auxiliary queue
insert source into the queue

while queue is not empty:
    node := front node in queue
    pop queue

    mark node as visited

    /*
        process the node
    */

    for each neighbor of node:
        if neighbor is not visited:
            insert neighbor to queue

```

<b>BFS on adjacency list</b>

```c++
queue<int> q;
q.push(source_node);

while(!q.empty()){
    int node = q.front();
    q.pop();

    visited[node] = true;

    for (int neighbor : adjList[node]){
        if (!visited[neighbor])
            q.push(neighbor);
    }
}
```

<b>BFS on adjacency matrix</b>

```c++
queue<int> q;
q.push(source_node);

while(!q.empty()){
    int node = q.front();
    q.pop();

    visited[node] = true;

    for (int neighbor=0;neightbor<N;neighbor++){
        if (!visited[neighbor] && adjMat[node][neighbor] != 0)
            q.push(neighbor);
    }
}
```

<b>BFS on grid</b>


```c++
//                                   Up     Down     Left     Right
vector<pair<int,int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

queue<pair<int,int>> q;
q.push(source_row, source_col);

while(!q.empty()){
    auto [curr_row, curr_col] = q.front();
    q.pop();

    visited[curr_row][curr_col] = true;

    for (auto [add_row, add_col] : directions){
        int next_row = add_row + curr_row;
        int next_col = add_col + curr_col;

        if (next_row < 0 || next_col < 0 || next_row >= ROW || next_col >= COL || visited[next_row][next_col])
            continue;

        q.push(next_row, next_col);
    }
}
```
### Application of Graph Traversal

- [Counting Rooms](https://cses.fi/problemset/task/1192)
- [Building Roads](https://cses.fi/problemset/task/1666)
- [Build Teams](https://cses.fi/problemset/task/1668)
- [Number of Islands](https://leetcode.com/problems/number-of-islands/)
- [Flight Routes Check](https://cses.fi/problemset/task/1682)
- [Fence Planning](http://www.usaco.org/index.php?page=viewproblem2&cpid=944)

### Using BFS to find shortest path

- [Message Route](https://cses.fi/problemset/task/1667)
- [Labyrinth](https://cses.fi/problemset/task/1193)
- [Monster](https://cses.fi/problemset/task/1194)
- [Network Delay time](https://leetcode.com/problems/network-delay-time/)
- [Construct a Palindrome](https://atcoder.jp/contests/abc197/tasks/abc197_f)

## Shortest Path Algorithms
### Dijkstra's Algorithm

<i>Pseudocode:</i>

```
distance[N] := [inf, inf, ...]
initialize an auxiliary Min-Heap

distance[source_node] := 0
insert a tuple of {0, source node} into the heap

while heap is not empty:
    node := top node in heap
    pop heap

    for each neighbor of node:
        if distance[neighbor] > distance[node] + weight between node and neighbor:
            distance[neighbor] := distance[node] + weight between node and neighbor
            insert a tuple of {distance[neighbor], neighbor node} into heap
```

<b>Dijkstra on adjacency list</b>

```c++
typedef pii pair<int,int>;
vector<int> distance(N, INT_MAX);
priority_queue<pii, vector<pii>, grater<pii>> pq;

distance[source_node] = 0;
pq.emplace(0, source_node);

while(!pq.empty()){
    auto [dist_to_node, node] = pq.top();
    pq.pop();

    for (auto [neighbor, weight] : adjList[node]){
        if (distance[neighbor] > distance[node] + weight){
            distance[neighbor] = distance[node] + weight;
            pq.emplace(distance[neighbor], neighbor);
        }
    }
}
```

<b>Dijkstra on adjacency matrix</b>

```c++
typedef pii pair<int,int>;
vector<int> distance(N, INT_MAX);
priority_queue<pii, vector<pii>, greater<pii>> pq;

distance[source_node] = 0;
pq.emplace(0, source_node);

while(!pq.empty()){
    auto [dist_to_node, node] = pq.top();
    pq.pop();

    for (int neighbor=0; neighbor<N; neighbor++){
        if (adjMat[node][neighbor] != -1 && distance[neighbor] > distance[node] + adjMat[node][neighbor]){
            distance[neighbor] = distance[node] + adjMat[node][neighbor];
            pq.emplace(distance[neighbor], neighbor);
        }
    }
}
```

<b>Dijkstra on grid</b>

```c++
//                                   Up     Down     Left     Right
vector<pair<int,int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

typedef tiii tuple<int, int, int>;
vector<vector<int>> distance(N, vector<int>(N, INT_MAX));
priority_queue<tiii, vector<tiii>, greater<tiii>> pq;

distance[source_row][source_col] = 0;
pq.emplace(0, source_row, source_col);

while(!pq.empty()){
    auto [dist_to_node, curr_row, curr_col] = pq.top();
    pq.pop();

    for (auto [add_row, add_col] : directions){
        int next_row = add_row + curr_row;
        int next_col = add_col + curr_col;

        if (next_row < 0 || next_col < 0 || next_row >= ROW || next_col >= COL)
            continue;

        if (distance[next_row][next_col] > distance[curr_row][curr_col] + grid[next_row][next_col]){
            distance[next_row][next_col] = distance[curr_row][curr_col] + grid[next_row][next_col]
            pq.emplace(distance[next_row][next_col], next_row, next_col);
        }
    }
}
```

### Bellman - Ford Algorithm

<i>Pseudocode:</i>

```
Bellman_Ford(source_node):
    distance[N] = [inf, inf, ...]
    distance[source_node] := 0

    Repeat N-1 times:
        for each edge in graph:
            get [u, v, weight] from edge
            if distance[v] > distance[u] + weight:
                distance[v] := distance[u] + weight
    
    for each edge in graph:
        get [u, v, weight] from edge
        if distance[v] > distance[u] + weight:
            log << "The graph has a negative cycle"
            return
```

<b>Bellman - Ford on edge list</b>

```c++
void Bellman_Ford(int source_node, vector<tuple<int,int,int>> &edgeList){
    vector<int> distance(N, 1e6);
    distance[source_node] = 0;

    for (int i=1;i<=N-1;i++){
        for (auto [u, v, weight] : edgeList){
            if (distance[v] > distance[u] + weight)
                distance[v] = distance[u] + weight; 
        }
    }

    for (auto [u, v, weight] : edgeList){
        if (distance[v] > distance[u] + weight){
            cout << "The graph has a negative cycle";
            return;
        }
    }
}
```

<b>Bellman - Ford on adjacency list</b>

```c++
void Bellman_Ford(int source_node, vector<pair<int,int>> &adjList){
    vector<int> distance(N, 1e6);
    distance[source_node] = 0;

    for (int i=1;i<=N-1;i++){
        for (int u=0;u<N;u++){
            for (auto [v, weight] : adjList[u]){
                if (distance[v] > distance[u] + weight)
                    distance[v] = distance[u] + weight; 
            }
        }
    }

    for (int u=0;u<N;u++){
        for (auto [v, weight] : adjList[u]){
            if (distance[v] > distance[u] + weight){
                cout << "The graph has a negative cycle";
                return;
            }
        }
    }
}
```

### Floyd - Warshall Algorithm

<i>Pseudocode:</i>

```
distance[N][N] := [[inf, inf, ...], [inf, inf, ...], ...]
for each edge in graph:
    get [u, v, weight] from edge
    distance[u][v] := weight

for each vertex v in graph:
    distance[v][v] := 0

for k from 0 to N-1:
    for u from 0 to N-1:
        for v from 0 to N-1:
            if distance[u][v] > distance[u][k] + distance[k][v]:
                distance[u][v] := distance[u][k] + distance[k][v]
```

<b>Floyd - Warshall on edge list</b>

```c++
vector<vector<int>> distance(N, vector<int>(N, 1e6));
for (auto [u, v, weight] : edgeList)
    distance[u][v] = weight;

for (int vertex=0;vertex<N;vertex++)
    distance[vertex][vertex] = 0;

for (int k=0;k<N;k++)
    for (int u=0;u<N;u++)
        for (int v=0;v<N;v++)
            if (distance[u][v] > distance[u][k] + distance[k][v])
                distance[u][v] = distance[u][k] + distance[k][v];
```

<b>Floyd - Warshall on adjacency list</b>

```c++
vector<vector<int>> distance(N, vector<int>(N, 1e6));
for (int u=0;u<N;u++)
    for (auto [v, weight] : adjList[u])
        distance[u][v] = weight;

for (int vertex=0;vertex<N;vertex++)
    distance[vertex][vertex] = 0;

for (int k=0;k<N;k++)
    for (int u=0;u<N;u++)
        for (int v=0;v<N;v++)
            if (distance[u][v] > distance[u][k] + distance[k][v])
                distance[u][v] = distance[u][k] + distance[k][v];
```

<b>Floyd - Warshall on adjacency matrix</b>

```c++
/*
    We use adjMat same as a distance array
    So any edges connect to the same vertex must be set to 0
    and any edges that do not exist must be set to inf
*/
for (int k=0;k<N;k++)
    for (int u=0;u<N;u++)
        for (int v=0;v<N;v++)
            if (adjMat[u][v] > adjMat[u][k] + adjMat[k][v])
```

### Application of Shortest Path Algorithms

- [Shortest Routes I](https://cses.fi/problemset/task/1671)
- [Shortest Routes II](https://cses.fi/problemset/task/1672)
- [Flight Discount](https://cses.fi/problemset/task/1195)
- [High Score](https://cses.fi/problemset/task/1673)
- [Cheapest Flights Within K Stops](https://leetcode.com/problems/cheapest-flights-within-k-stops)
- [Cycle Finding](https://cses.fi/problemset/task/1197)
- [Longest Flight Route](https://cses.fi/problemset/task/1680)

## Minimum Spanning Tree Algorithms

### Union Find (Disjoint Set) structure

```c++
struct UnionFind {
    //Assume that a node is represented as an interger
    int rank[N];
    int parent[N];

    UnionFind(){
        for (int i=0;i<N;i++){
            rank[i] = 0;
            parent[i] = i;
        }
    }

    int find(int u){
        if (u == parent[u]) return parent[u];
        return parent[u] = find(parent[u]);
    }

    void union(int u, int v){
        int pu = parent[u], pv = parent[v];
        if (rank[pu] < rank[pv]) parent[pu] = pv;
        else if (rank[pu] > rank[pv]) parent[pv] = pu;
        else{
            parent[pu] = pv;
            rank[pv]++;
        }
    }
};
```

### Kruskal's Algorithm

<i>Pseudocode:</i>

```
sort all edges according to their weight
initialize a set with Union Find structure

for each edge in graph:
    get [u, v, weight] from edge
    if u and v are not in the same set:
        union set of u and v
        add the current edge to the tree
```

<b>Kruskal on edge list</b>

```c++
sort(edgeList.begin(), edgeList.end());
UnionFind UF;
for (auto &[u, v, weight] : edgeList){
    if (UF.find(u) != UF.find(v)){
        UF.union(u, v);
        /*Add the current edge to the tree*/
    }
}
```

### Prim's Algorithm

<i>We don't talk about Prim's around here boi</i> 🤫

### Application of Minimum Spanning Tree

- [Road Reparation](https://cses.fi/problemset/task/1675)
- [Superbull](http://www.usaco.org/index.php?page=viewproblem2&cpid=531)
- [Min Cost to Connect All Points](https://leetcode.com/problems/min-cost-to-connect-all-points/)
- [กระเช้าไฟฟ้า](https://programming.in.th/tasks/toi12_cablecar)