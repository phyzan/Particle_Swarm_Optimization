#include <condition_variable>
#include <filesystem>
#include <fstream>

#include "grid_search.hpp"

void init_folder(std::string name)
{
    std::filesystem::path folder(name);

    if (std::filesystem::is_directory(folder) == true)
    { // If the folder named exists, delete it and its contents.
        std::filesystem::remove_all(folder);
    }

    if (std::filesystem::create_directories(folder))
    { // Create the folder.
        std::cout << "Folder created successfully.\n";
    }
    else
    {
        std::cerr << "Failed to create folder.\n";
    }
}

std::string leading_zeros(int grid_size, int current)
{
    // Used to align the numbers of exited threads.

    int grid_length = int(std::to_string(grid_size).length());
    int current_length = int(std::to_string(current).length());
    std::string zero = "0";
    std::string res = "";

    if (current_length < grid_length)
    {
        for (int i = 0; i < (grid_length - current_length); i++)
        {
            res += "0";
        }
    }

    return res;
}

MPFR_ARR *create_grid_bounds(MPFR_ARR original_bounds, int res, int &grid_size)
{
    int dim = original_bounds.rows();                // Dimensions of the hyperspace.
    grid_size = pow(2, dim * res);                   // The number of subspaces that will be created.
                                                     // If the dimensions of the space are 2 and the
                                                     // resolution 3, 2^(2*3) = 64 new subspaces will
                                                     // be created.
    int dim_divisions = pow(2, res);                 // The number of times each dimension will be
                                                     // divided. For example if dim = 2 and res = 3,
                                                     // the x and y axis will be divided in 2^3 = 8
                                                     // each, in order to create an 8x8 grid, with
                                                     // grid_size = 2^(2*3) = 8*8 = 64 subspaces.
    MPFR_ARR dim_bounds(dim, dim_divisions + 1);     // An array to store the values of the new boundaries.
    MPFR_ARR *grid_bounds = new MPFR_ARR[grid_size]; // An array to store the new subspace boundaries.

    for (int i = 0; i < grid_size; i++)
    { // Initialize the arrays.
        grid_bounds[i] = MPFR_ARR(dim, 2);
    }

    for (int i = 0; i < dim; i++)
    { // Calculate and store the subspace size for each dimension.
        MP_REAL step = (original_bounds(i, 1) - original_bounds(i, 0)) / dim_divisions;

        for (int j = 0; j <= dim_divisions; j++)
        {
            dim_bounds(i, j) = original_bounds(i, 0) + step * j;
        }
    }

    for (int i = 0; i < dim; i++)
    { // O(dim * grid_size).
        int cell = 0;

        for (int j = 0; j < pow(dim_divisions, dim - i - 1); j++)
        { // This nested loop is not O(n^3).
            for (int bound = 0; bound < dim_divisions; bound++)
            { // It will always be O(grid_size).
                for (int k = 0; k < (grid_size / pow(dim_divisions, dim - i)); k++)
                {                                                       // The three loops are needed to repeat
                    grid_bounds[cell](i, 0) = dim_bounds(i, bound);     // each dimension's boundaries in different
                    grid_bounds[cell](i, 1) = dim_bounds(i, bound + 1); // intervals, to create and store all of
                                                                        // the subspaces' boundaries. It has to be
                    cell++;                                             // done in this way in order for the
                } // method to be scalable in n-dimensions.
            }
        }
    }

    return grid_bounds;
}

MPFR_ARR grid_search(grid_params gp, double swap_point, int precision, std::string Obj_F, bool constriction)
{
    int grid_size = 0;
    int running_tasks = 0;
    int launched_tasks = gp.resume_from_subspace;
    MPFR_ARR grid_minima; // An array that stores all the minima found.
    MPFR_EMPTY empty;
    std::mutex mutex;

    std::streambuf *prompt = std::cout.rdbuf();
    std::condition_variable exit_signal;
    std::vector<std::thread> threads;

    std::cout << std::setprecision(16);

    if (gp.type != "Global_Classic" && gp.type != "Local_Classic")
    {
        std::cout << "~> Error: Please select a correct type. [Global_Classic] or [Local_Classic]." << std::endl;

        return empty;
    }

    grid_minima = empty;

    std::cout << "Creating Grid Space... ";

    MPFR_ARR *grid_bounds = create_grid_bounds(gp.p.bounds, gp.resolution, grid_size);

    std::cout << "Complete." << std::endl;
    std::cout << "Running " << grid_size << " threads, in batches of " << gp.batch_size << ".\n" << std::endl;

    if (gp.resume_from_subspace > 0 && grid_size > gp.resume_from_subspace)
    {
        std::cout << "Resuming from subspace... " << gp.resume_from_subspace << std::endl;
    }
    else if (gp.resume_from_subspace > 0 && grid_size <= gp.resume_from_subspace)
    {
        std::cout << "Cannot resume from subspace... " << gp.resume_from_subspace << std::endl;
        std::cout << "Subspace number " << gp.resume_from_subspace << " is larger than grid size " << grid_size << "."
                  << std::endl;

        delete[] grid_bounds;
        grid_bounds = nullptr;

        return MPFR_EMPTY();
    }

    init_folder(gp.dest);

    auto thread_task = [&](int id, pso_params<MP_REAL, MPFR_ARR> p, enhanced_params<MP_REAL> ep,
                           obj_params<MP_REAL> objp, int precision, std::string type) {
        MPFR_ARR min; // The minimum found in this thread.
        std::string file_name = gp.dest + "/Grid_Search_Global_Classic_PSO_" + std::to_string(id) +
                                ".txt"; // The output file of this thread.
        std::ofstream output(file_name);

        if (type == "Global_Classic")
        {
            Global_Classic *global_classic_pso = Global_Classic_Init(p, swap_point, precision, &output);

            if (global_classic_pso != nullptr)
            {                                       // If the method is initialized correctly,
                global_classic_pso->print_params(); // print the parameters,

                global_classic_pso->set_Obj_F(Obj_F, objp); // set the objective function

                min = global_classic_pso->fit(); // save the minimum and

                mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);

                delete global_classic_pso; // delete the method instance.
                global_classic_pso = nullptr;
            }
        }
        else if (type == "Local_Classic")
        {
            Local_Classic *local_classic_pso = Local_Classic_Init(p, ep, swap_point, precision, constriction, &output);

            if (local_classic_pso != nullptr)
            {
                local_classic_pso->print_params();

                local_classic_pso->set_Obj_F(Obj_F, objp);

                min = local_classic_pso->fit();

                mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);

                delete local_classic_pso;
                local_classic_pso = nullptr;
            }
        }

        {
            std::lock_guard<std::mutex> lock(mutex); // One thread allowed to enter at a time.

            std::cout.rdbuf(prompt); // Redirect to terminal output.

            if (grid_minima.size() == 0 && min.size() != 0)
            {                      // If grid_minima is empty,
                grid_minima = min; // just save the found minimum.
            }
            else if (min.size() != 0)
            {                                                                     // else,
                MPFR_ARR new_grid_minima(p.dim, grid_minima.cols() + min.cols()); // expand the grid_minima and

                new_grid_minima.leftCols(grid_minima.cols()) = grid_minima; // save the new found minimum.
                new_grid_minima.rightCols(min.cols()) = min;

                grid_minima = new_grid_minima;
            }

            std::cout << "Thread " << leading_zeros(grid_size, id) << id << " exited." << std::endl;

            running_tasks--;
        }

        exit_signal.notify_one();
    };

    while (launched_tasks < grid_size)
    { // Loop that keeps batch_size number of tasks running.
        std::unique_lock<std::mutex> lock(mutex);

        exit_signal.wait(lock, [&] { return running_tasks < gp.batch_size; });

        gp.p.bounds = grid_bounds[launched_tasks];

        threads.emplace_back(thread_task, launched_tasks, gp.p, gp.ep, gp.objp, precision, gp.type);

        launched_tasks++;
        running_tasks++;
    }

    for (auto &t : threads)
    { // Wait for all threads to finish.
        t.join();
    }

    delete[] grid_bounds;
    grid_bounds = nullptr;

    return grid_minima;
}
