#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <memory>
#include <optional>
#include <random>
#include <vector>

// Random device seed
constexpr auto kSeed = 1234;

constexpr float kEpsilon = 1e-6;

constexpr float kXmin = -1.0;
constexpr float kXmax = 1.0;
constexpr float kYmin = -1.0;
constexpr float kYmax = 1.0;

constexpr float kMinDistanceX = 1e-5;
constexpr float kMinDistanceY = 1e-5;

constexpr int kPointsNumber = 10;  // 0000;

class TimeBench {
 public:
  TimeBench(std::string aLabel) : mLabel{std::move(aLabel)} {}

  ~TimeBench() {
    const auto endPoint = std::chrono::high_resolution_clock::now();
    const auto processTime =
        std::chrono::duration_cast<std::chrono::milliseconds>(endPoint -
                                                              mStartPoint)
            .count();
    std::cout << mLabel << " : " << processTime << " milliseconds\n";
  }

 private:
  const std::string mLabel;
  const std::chrono::high_resolution_clock::time_point mStartPoint =
      std::chrono::high_resolution_clock::now();
};

struct Point {
  Point(float aX, float aY) : x{aX}, y{aY} {}
  float x{};
  float y{};

  float operator[](uint aCurrentDimension) const {
    switch (aCurrentDimension % 2) {
      case 0:
        return x;
        break;
      case 1:
        return y;
        break;
    }
    return NAN;
  }
};

bool operator==(const Point& lhs, const Point& rhs) {
  return (std::abs(lhs.x - rhs.x) < kEpsilon &&
          std::abs(lhs.y - rhs.y) < kEpsilon);
}

bool operator!=(const Point& lhs, const Point& rhs) { return !(lhs == rhs); }

std::ostream& operator<<(std::ostream& os, const Point& data) {
  os << "{" << data.x << "," << data.y << "}";
  return os;
}

struct Rectangle {
  Rectangle(float aXmin, float aXmax, float aYmin, float aYmax)
      : mXmin{std::min(aXmin, aXmax)},
        mXmax{std::max(aXmin, aXmax)},
        mYmin{std::min(aYmin, aYmax)},
        mYmax{std::max(aYmin, aYmax)} {}

  [[nodiscard]] bool isPointInside(const Point& aPoint) const {
    if (aPoint.x > mXmax || aPoint.x < mXmin || aPoint.y > mYmax ||
        aPoint.y < mYmin) {
      // std::cout << "Out of borders" << '\n';
      return false;
    }
    return true;
  }

  operator bool() const {
    return !(std::isfinite(mXmin) && std::isfinite(mXmax) &&
             std::isfinite(mYmin) && std::isfinite(mYmax));
  }

  [[nodiscard]] float Xmin() const { return mXmin; }

  [[nodiscard]] float Xmax() const { return mXmax; }

  [[nodiscard]] float Ymin() const { return mYmin; }

  [[nodiscard]] float Ymax() const { return mYmax; }

  float Max(uint aCurrentDimension) const {
    switch (aCurrentDimension % 2) {
      case 0:
        return mXmax;
        break;
      case 1:
        return mYmax;
        break;
    }
    return NAN;
  }

  float Min(uint aCurrentDimension) const {
    switch (aCurrentDimension % 2) {
      case 0:
        return mXmin;
        break;
      case 1:
        return mYmin;
        break;
    }
    return NAN;
  }

 private:
  // Borders
  float mXmax;
  float mXmin;
  float mYmax;
  float mYmin;
};

std::ostream& operator<<(std::ostream& os, const Rectangle& data) {
  os << "{" << data.Xmin() << "," << data.Xmax() << "}{" << data.Ymin() << ","
     << data.Ymax() << "}";
  return os;
}

class Node {
 public:
  static std::unique_ptr<Node> createNode(const Point& aPoint) {
    std::unique_ptr<Node> result(new Node(aPoint));
    return result;
  }

  bool insertPoint(const Point& aPoint, uint aInsertionDepth = 0) {
    std::cout << "Insert " << aPoint << " at depth " << aInsertionDepth << '\n';

    // Check if point is previously inserted
    if (aPoint == mPoint) {
      if (mIsdeleted) {
        std::cout << "Restore deleted point\n";
        mIsdeleted = false;
        return true;
      }
      std::cout << "Point already exists!\n";
      return false;
    }

    bool isInserted = true;

    uint aCurrentDimension = aInsertionDepth % 2;

    if (aPoint[aCurrentDimension] < mPoint[aCurrentDimension]) {
      std::cout << "Insert Left\n";
      if (mLeft) {
        isInserted = mLeft->insertPoint(aPoint, aInsertionDepth + 1);
      } else {
        mLeft = createNode(aPoint);
        isInserted = static_cast<bool>(mLeft);
      }
    } else {
      std::cout << "Insert Right\n";
      if (mRight) {
        isInserted = mRight->insertPoint(aPoint, aInsertionDepth + 1);
      } else {
        mRight = createNode(aPoint);
        isInserted = static_cast<bool>(mRight);
      }
    }

    if (isInserted) {
      updateDepth();
    }

    return isInserted;
  }

  bool findPoint(const Point& aPoint, uint aSearchDepth = 0) {
    std::cout << "Find " << aPoint << " at depth " << aSearchDepth << '\n';

    // Check if Node point found is close enough
    if (mPoint == aPoint && !mIsdeleted) {
      std::cout << "Found\n";
      return true;
    }

    // Search children
    uint aCurrentDimension = aSearchDepth % 2;

    if (aPoint[aCurrentDimension] < mPoint[aCurrentDimension]) {
      std::cout << "Find Left\n";
      if (mLeft) {
        return mLeft->findPoint(aPoint, aSearchDepth + 1);
      }
    } else {
      std::cout << "Find Right\n";
      if (mRight) {
        return mRight->findPoint(aPoint, aSearchDepth + 1);
      }
    }

    std::cout << "Not Found\n";
    return false;
  }

  [[nodiscard]] std::list<Point> findPointsInArea(const Rectangle& aArea,
                                                  uint aSearchDepth = 0) const {
    std::cout << "Find in " << aArea << " at depth " << aSearchDepth << '\n';

    std::list<Point> result;

    if (!mIsdeleted && aArea.isPointInside(mPoint)) {
      result.push_back(mPoint);
    }

    // Search children
    // TODO: Parallelization point
    uint aCurrentDimension = aSearchDepth % 2;

    if (mPoint[aCurrentDimension] > aArea.Min(aCurrentDimension) && mLeft) {
      std::cout << "Search Area Left\n";
      auto&& pointsFoundInArea =
          mLeft->findPointsInArea(aArea, aSearchDepth + 1);
      result.splice(std::end(result), pointsFoundInArea);
    }

    if (mPoint[aCurrentDimension] < aArea.Max(aCurrentDimension) && mRight) {
      std::cout << "Search Area Right\n";
      auto&& pointsFoundInArea =
          mRight->findPointsInArea(aArea, aSearchDepth + 1);
      result.splice(std::end(result), pointsFoundInArea);
    }

    return result;
  }

  bool deletePoint(const Point& aPoint, uint aDeleteDepth = 0) {
    std::cout << "Delete " << aPoint << " at depth " << aDeleteDepth << '\n';

    // Check if Node point to delete is close enough
    if (!mIsdeleted && mPoint == aPoint) {
      mIsdeleted = true;
      // Current node point value is invalid - update depth
      const auto largestDepth = std::max({getDepth(mLeft), getDepth(mRight)});
      mDepth = largestDepth;

      std::cout << (mDepth ? "Deleted current node"
                           : "Current node invalidated")
                << '\n';

      return true;
    }

    // Delete in children
    uint aCurrentDimension = aDeleteDepth % 2;

    if (aPoint[aCurrentDimension] < mPoint[aCurrentDimension]) {
      std::cout << "Delete Left\n";
      if (mLeft) {
        const auto isDeleted = mLeft->deletePoint(aPoint, aDeleteDepth + 1);
        if (mLeft->isEmpty()) {
          mLeft.reset();
        }
        if (isDeleted) {
          updateDepth();
        }
        return isDeleted;
      }
    } else {
      std::cout << "Delete Right\n";
      if (mRight) {
        const auto isDeleted = mRight->deletePoint(aPoint, aDeleteDepth + 1);
        if (mRight->isEmpty()) {
          mRight.reset();
        }
        if (isDeleted) {
          updateDepth();
        }
        return isDeleted;
      }
    }

    std::cout << "Point " << aPoint << " not deleted\n";
    return false;
  }

  [[nodiscard]] Point getPoint() const { return mPoint; }

  [[nodiscard]] std::list<Point> getAllPoints(uint aSearchDepth = 0) const {
    std::list<Point> points;

    if (!mIsdeleted) {
      points.push_back(mPoint);
    }

    // Search children
    uint aCurrentDimension = aSearchDepth % 2;

    // TODO: Parallelization point
    if (mLeft) {
      std::cout << "Search Area Left\n";
      auto&& pointsFound = mLeft->getAllPoints(aSearchDepth + 1);
      points.splice(std::end(points), pointsFound);
    }

    if (mRight) {
      std::cout << "Search Area Right\n";
      auto&& pointsFound = mRight->getAllPoints(aSearchDepth + 1);
      points.splice(std::end(points), pointsFound);
    }

    return points;
  }

  [[nodiscard]] int getDepth() const { return mDepth; }

  [[nodiscard]] bool isEmpty() const { return mDepth == 0; }

 private:
  Node(const Point& aPoint) : mPoint{aPoint} {
    // std::cout << "Node point " << mPoint << '\n';
  }

  static int getDepth(std::unique_ptr<Node>& aNode) {
    if (aNode) {
      return aNode->mDepth;
    }
    return 0;
  }

  void updateDepth() {
    mDepth = std::max({getDepth(mLeft), getDepth(mRight)});
    if (!mIsdeleted) {
      mDepth++;
    }
  }

  // Point
  Point mPoint;

  // Node depth
  int mDepth{1};

  // TODO: Refactor using mDepth
  bool mIsdeleted{false};

  // Node quads
  std::unique_ptr<Node> mLeft;
  std::unique_ptr<Node> mRight;
};

int main(int /*argc*/, char* /*argv*/[]) {
  std::vector<Point> testPoints;
  testPoints.reserve(kPointsNumber);

  // Random generator of points
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(kSeed);

  std::uniform_real_distribution<float> distX{kXmin, kXmax};
  std::uniform_real_distribution<float> distY{kYmin, kYmax};

  std::cout << "Points count: " << kPointsNumber << '\n';

  // Insert points in vector
  {
    TimeBench bench{"Insert points in vector"};
    for (int i = 0; i < kPointsNumber; i++) {
      testPoints.emplace_back(distX(gen), distY(gen));
    }
  }

  // Root node
  auto root = Node::createNode(testPoints.at(0));
  if (!root) {
    std::cout << "Failed to create root!\n";
    return EXIT_FAILURE;
  }

  // Insert points in KDTree
  int insertedPoints = 1;
  {
    TimeBench bench{"Insert points in KDTree"};

    // First point already inserted
    for (auto startIt = std::next(std::begin(testPoints));
         startIt != std::end(testPoints); ++startIt) {
      const auto point = *startIt;
      // std::cout << "===\n";
      const auto isInserted = root->insertPoint(point);
      if (isInserted) {
        insertedPoints++;
        // std::cout << "Root depth after insert: " << root->getDepth() << '\n';
      } else {
        // std::cout << "Failed to insert point!" << '\n';
      }
      // std::cout << "===\n";
    }
  }
  if (insertedPoints != testPoints.size()) {
    std::cout << "Expected points: " << testPoints.size()
              << "  Inserted points: " << insertedPoints << '\n';
  }

  // Find points in vector
  int pointsFound = 0;
  {
    TimeBench bench{"Find points in vector"};
    for (const auto& point : testPoints) {
      auto itFoundPoint =
          std::find(std::begin(testPoints), std::end(testPoints), point);
      if (std::end(testPoints) != itFoundPoint) {
        pointsFound++;
      }
    }
  }
  if (pointsFound != testPoints.size()) {
    std::cout << "Expected points: " << testPoints.size()
              << "  Found points: " << pointsFound << '\n';
  }

  // Find points in KDTree
  pointsFound = 0;
  {
    TimeBench bench{"Find points in KDTree"};
    for (const auto& point : testPoints) {
      std::cout << "=====\n";
      const auto isPointFound = root->findPoint(point);
      std::cout << "Find point: " << isPointFound << '\n';
      if (isPointFound) {
        pointsFound++;
      }
      std::cout << "=====\n";
    }
  }
  if (pointsFound != testPoints.size()) {
    std::cout << "Expected points: " << testPoints.size()
              << "  Found points: " << pointsFound << '\n';
  }

  // Missing point?
  // const auto pointFound = root->findPoint(Point{0.1, 0.1});
  // std::cout << "Find point: " << pointFound << '\n';

  // Search points in area
  const Rectangle searchArea{-1.0, 0.0, -1.0, 0.0};
  std::vector<Point> pointsFoundInArea;

  // Search points in area in vector
  {
    TimeBench bench{"Find points in area in vector"};
    std::copy_if(std::begin(testPoints), std::end(testPoints),
                 std::back_inserter(pointsFoundInArea),
                 [&searchArea](const auto& aPoint) {
                   return searchArea.isPointInside(aPoint);
                 });
  }
  std::cout << "Points found in area: " << pointsFoundInArea.size() << '\n';

  // Search points in area in KDTree
  std::list<Point> pointsFoundInAreaQT;
  {
    TimeBench bench{"Find points in area in KDTree"};
    pointsFoundInAreaQT = root->findPointsInArea(searchArea);
  }
  std::cout << "Points found in area: " << pointsFoundInAreaQT.size() << '\n';

  // Delete points in KDTree
  int pointsDeleted = 0;
  {
    TimeBench bench{"Delete points in KDTree"};
    for (const auto& point : testPoints) {
      std::cout << "=====\n";
      std::cout << "Points in root:\n";
      for (const auto& currentPoint : root->getAllPoints()) {
        std::cout << currentPoint << '\n';
      }
      std::cout << "=====\n";
      const auto isPointDeleted = root->deletePoint(point);
      if (isPointDeleted) {
        pointsDeleted++;
      }
      std::cout << "Delete point: " << point << " " << isPointDeleted
                << "  root depth: " << root->getDepth() << '\n';
      std::cout << "=====\n";
    }
  }

  if (pointsDeleted != testPoints.size()) {
    std::cout << "Expected deleted points: " << testPoints.size()
              << "  Deleted points: " << pointsDeleted << '\n';
  }

  std::cout << "Root empty: " << root->isEmpty() << '\n';

  std::cout << "=====\n";
  std::cout << "Points in root:\n";
  for (const auto& currentPoint : root->getAllPoints()) {
    std::cout << currentPoint << '\n';
  }
  std::cout << "=====\n";

  std::cout << "Done.\n";

  return EXIT_SUCCESS;
}