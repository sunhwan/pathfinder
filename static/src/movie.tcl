# Generates movie from pathway
# This script is intended for VMD.

mol load pdb pathway.pdb
mol modstyle 0 top tube

        # Loop through every frame in the trajectory.
        #
        for { set frame 0 } { $frame < [ molinfo top get numframes ] } { incr frame } {

                # Set the current frame.
                #
                animate goto $frame

                # Render snapshot of the system at the current frame.
                #
                set filename frame.[ format "%05d" $frame ].rgb
                render TachyonInternal $filename
        }

exit
