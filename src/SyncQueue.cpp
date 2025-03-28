#include "SyncQueue.h"

void SyncQueue::push(const std::pair<std::vector<dy4::real>, std::vector<dy4::real>> &data) {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [this] { return queue.size() < maxSize || done; }); // Wait if the queue is full

    if (done) return; // Exit if the queue is marked as done

    queue.push(data);
    cv.notify_one(); // Notify a waiting thread
}

bool SyncQueue::pop(std::pair<std::vector<dy4::real>, std::vector<dy4::real>> &data) {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [this] { return !queue.empty() || done; });

    if (queue.empty()) {
        return false;
    }

    data = queue.front();
    queue.pop();
    cv.notify_one();
    return true;
}

void SyncQueue::setDone() {
    {
        std::lock_guard<std::mutex> lock(mtx);
        done = true;
    }
    cv.notify_all();
}

bool SyncQueue::isDone() {
    std::lock_guard<std::mutex> lock(mtx);
    return done;
}