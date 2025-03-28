#ifndef SYNC_QUEUE_H
#define SYNC_QUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>
#include <utility>
#include <vector>
#include "dy4.h"

// using namespace std;

class SyncQueue {
public:
    void push(const std::pair<std::vector<dy4::real>, std::vector<dy4::real>> &data);
    bool pop(std::pair<std::vector<dy4::real>, std::vector<dy4::real>> &data);
    void setDone();
    bool isDone();

private:
    std::queue<std::pair<std::vector<dy4::real>, std::vector<dy4::real>>> queue;
    std::mutex mtx;
    std::condition_variable cv;
    bool done = false;
    size_t maxSize = 5; // Maximum size of the queue
};

#endif // SYNC_QUEUE_H
