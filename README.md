# University of Toronto's PHY254 "Classical Mechanics" lecture notes, summer 2023 edition

## On authorship

Lecture notes are not like articles or books, but they are still protected by copyright. I cannot claim sole authorship on these lecture notes. Many times, I drew inspiration from previous colleagues who taught this course. This is usually fine, as lecture notes are passed on in person and remain hidden from public view.

Here however, to ease my student's access to these notes, I decided to make this repository public, while sometimes re-using the exact same words of my previous colleagues. I therefore need to acknoledge my "co-authors":
* [Paul Kushner](https://www.pjk.atmosp.physics.utoronto.ca/) (taught PHY254 in 2009 and 2011),
* [Sabine Stanley](https://sabinestanley.com/), now at Johns Hopkins University (taught PHY254 in 2010 and 2012),
* [Stephen Morris](https://imgflip.com/memegenerator/162372564/Domino-Effect), happily retired (taught PHY254 in 2013-2015 and 2020).

If you recognize some of your prose in these notes, please reach out to be added on the list.

## Using this repository with the University of Toronto's JupyterHub

Note that most of these instructions would work more or less the same way locally on Linux and MacOS (I have no idea about Windows) and advanced users will not need my instructions, nor our University's JupyterHub, to use these lectures notes. For all others, here is a solution that should work without needing to configure your computer, as long as you have internet access. Please contact me if these instructions don't work, or click on the "Issues" tab above.

### If first cloning into the JupyterHub:

1. Go to https://jupyter.utoronto.ca, log in with your UTorID, and start the server.
2. Near the top-right-hand corner of the home menu, hit the drop-down menu "New", and click on Terminal.
3. Clone the repo of my chapters. There are [two possibilities](https://www.howtogeek.com/devops/should-you-use-https-or-ssh-for-git/):
    * The **SSH protocol** is more secure but requires you to first [generate a new SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key), applying these instructions in your JupyterHub terminal. Once this is sorted, enter the command
        ```
        git clone git@github.com:ngrisouard/254-lectures-2023S.git
        ```
        in the Terminal and hit return. 

    * An alternative is to use the **HTTP protocol**.
        ```
        git clone https://github.com/ngrisouard/254-lectures-2023S.git
        ```
        To be honest, I don't use this protocol and I don't know if the rest of my instructions still apply if you use this protocol.
      
   Either options should create a new directory called `254-lectures-2023S`, containing my lecture notes.

### Refreshing the lecture notes

At the beginning of the term, the lectures will be quite empty. As the term progresses, I will add to them, correct typos, etc. You will want to refresh the contents periodically, before and after each lecture.

4. To get the latest updates, repeat steps 1-2 above.
5. In the Terminal, change directory to go to the repository by typing
    ```
    cd 254-lectures-2023S
    ```
    and hit return.

6. If you just want to play around with the notes a little in-between two updates but want your repository to closely match what I have on GitHub, enter the commands
    ```
    git fetch --all
    git reset --hard origin/master
    ```
    in the Terminal, where you left at the end of step 5, whenever you want to re-align your notes with mine.

7. If you want to be more fancy, modify the notebooks significantly and keep your modifications while keeping up-to-date with my content, I suggest you create your own "fork" and work on it. You can then sync your fork with my repository by following these instructions: https://help.github.com/en/articles/syncing-a-fork

### Opening a Jupyter notebook

8. Repeat steps 1-3 above.
9. Repeat step 6 or 7 above, probably.
10. You can navigate to the Jupyter file, using the graphical interface of the home menu of syzygy. You are looking for a `.ipynb` file.

**Please let me know if this procedure does not work.**
