import tkinter as tk
from tkinter import ttk, filedialog, Listbox, Text, END, messagebox
import threading
import queue # For communication between threads
import os # For path manipulation (basename, splitext)
import re # For regular expressions to extract sort keys
import ram_api_v1 as ram_api_module # Make sure this import matches your file/module name

# --- Helper function for custom sorting ---
def extract_sort_parts_for_file_list(filepath):
    """
    Extracts sorting keys (letter, number) from a filepath for custom sorting.
    Sorts by:
    1. A representative letter (descending: Z before A).
       - Derived from the last character of the last word in the filename stem.
    2. A representative number (descending: 10 before 1).
       - Derived from the last sequence of digits in the filename stem.

    Args:
        filepath (str): The full path to the file.

    Returns:
        tuple: (letter_key, num_key)
               letter_key is a character (e.g., 'Z', 'A', or '' if no letters).
               num_key is an integer (e.g., 10, 1, or -1 if no numbers).
    """
    try:
        name_part = os.path.splitext(os.path.basename(filepath))[0]

        # Extract number key
        numbers_found = re.findall(r'\d+', name_part)
        num_key = int(numbers_found[-1]) if numbers_found else -1 # Default for no numbers

        # Extract letter key
        words_found = re.findall(r'[A-Za-z]+', name_part)
        letter_key = '' # Default for no letters
        if words_found:
            last_word = words_found[-1]
            if last_word: # Ensure last_word is not empty
                letter_key = last_word[-1].upper() # Last char of last word, uppercased

        return (letter_key, num_key)
    except Exception:
        # Fallback for any parsing error, ensures sorting doesn't crash
        # Files causing errors will likely sort together at one end.
        return ('', -1)


# --- Actual RAM Concept Processing Logic ---
# This function will be run in a separate thread.
# It now assumes that ram_api functions accept 'progress_queue' as an argument.
def ram_load_rundown_processing(file_path_list, progress_queue):
    """
    Processes the RAM Concept files using the provided ram_api.
    This function is run in a separate thread.
    It assumes ram_api.delete_exisitng_loads and ram_api.ram_load_rundown
    have been modified to accept progress_queue as an argument for detailed feedback.

    Args:
        file_path_list (list): A list of file paths in the desired order.
        progress_queue (queue.Queue): A queue to send progress messages to the UI.
    """

    try:
        # Append path and import ram_api within the thread
        # IMPORTANT: Ensure this path is correct for your system
        if not file_path_list or len(file_path_list) < 2:
            progress_queue.put("Error: At least two files are required for the rundown process.\n")
            progress_queue.put("COMPLETED_WITH_ERROR")
            return

        progress_queue.put("Starting RAM Concept load rundown...\n")
        if hasattr(ram_api_module, '__file__') and ram_api_module.__file__: # Check if __file__ exists and is not None
            progress_queue.put("Opening RAM Concept... \n")
        else:
            progress_queue.put("Opening RAM Concept (unable to determine file location).\n")

        ram_concept_wrapper_instance = ram_api_module.RamConcept()
        ram_concept_wrapper_instance.create_concept()
        concept = ram_concept_wrapper_instance.concept

        for i in range(len(file_path_list) - 1):
            reaction_path = file_path_list[i]
            target_path = file_path_list[i + 1]

            progress_queue.put(f"\n--- Processing Pair {i+1} of {len(file_path_list)-1} ---\n")
            progress_queue.put(f"Reaction file (source): {reaction_path}\n")
            progress_queue.put(f"Target file (destination): {target_path}\n")
            
            # Call delete_exisitng_loads, passing the progress_queue
            ram_api_module.delete_existing_loads(target_path, concept, progress_queue)

            # Call ram_load_rundown, passing the progress_queue
            ram_api_module.ram_load_rundown(reaction_path, target_path, concept, progress_queue)
            
            progress_queue.put(f"--- Pair {i+1} completed. ---\n")

        progress_queue.put("\nShutting down RAM Concept...\n")
        ram_concept_wrapper_instance.shutdown_concept()
        
        progress_queue.put("\nLoad rundown process completed successfully for all pairs!\n")
        progress_queue.put("ALL_DONE")

    except ImportError as ie:
        progress_queue.put(f"ImportError: Could not import 'ram_api'.\n")
        progress_queue.put(f"Please ensure 'ram_concept_wrapper_instance_tools' is in the specified path and the module/file name is correct.\n")
        progress_queue.put(f"Current sys.path includes: C:\\_Github\\structural_engineering_toolbox\\ram_concept_wrapper_instance_tools\n")
        progress_queue.put(f"Details: {str(ie)}\n")
        progress_queue.put("COMPLETED_WITH_ERROR")
    except (AttributeError, TypeError) as ae_te: 
        progress_queue.put(f"Error calling API function ({type(ae_te).__name__}):\n")
        progress_queue.put(f"This can occur if a function is missing from 'ram_api', called incorrectly, \n"
                           f"or if 'ram_load_rundown' or 'delete_exisitng_loads' were not updated \n"
                           f"to accept 'progress_queue' as an argument.\n")
        progress_queue.put(f"Details: {str(ae_te)}\n")
        progress_queue.put("COMPLETED_WITH_ERROR")
    except Exception as e:
        progress_queue.put(f"An unexpected error occurred during RAM Concept processing:\n")
        progress_queue.put(f"Error type: {type(e).__name__}\n")
        progress_queue.put(f"Details: {str(e)}\n")
        progress_queue.put("COMPLETED_WITH_ERROR")
    finally:
        # This block will ALWAYS execute, whether an exception occurred or not.
        # This is where shutdown should happen.
        if ram_concept_wrapper_instance:
            progress_queue.put("\nEnsuring RAM Concept session is shut down...\n")
            # Call the shutdown method of your RamConcept wrapper class
            ram_concept_wrapper_instance.shutdown_concept() 
            progress_queue.put("RAM Concept shutdown process completed by wrapper.\n")
        else:
            progress_queue.put("\nRAM Concept wrapper was not initialized; no shutdown performed by UI thread.\n")
        

class LoadRundownApp:
    def __init__(self, master):
        self.master = master
        master.title("CED RAM Concept Load Rundown")
        master.geometry("1280x720") 

        style = ttk.Style()
        style.configure("TButton", padding=6, relief="flat", font=('Helvetica', 10))
        style.configure("TLabel", font=('Helvetica', 10))
        style.configure("Header.TLabel", font=('Helvetica', 12, 'bold'))

        self.selected_files = []
        self.progress_queue = queue.Queue() 
        self.processing_thread = None

        self.create_widgets()
        self.master.protocol("WM_DELETE_WINDOW", self.on_closing) 

    def create_widgets(self):
        # --- Main PanedWindow for resizable sections ---
        main_paned_window = ttk.PanedWindow(self.master, orient=tk.VERTICAL)
        main_paned_window.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # --- Top Frame: File Selection and Reordering ---
        top_frame = ttk.Frame(main_paned_window, padding=10)
        main_paned_window.add(top_frame, weight=1) 

        self.select_button = ttk.Button(top_frame, text="Select RAM Concept Files (.cpt)", command=self.select_files_action)
        self.select_button.pack(pady=(0, 10), fill=tk.X)

        list_controls_frame = ttk.Frame(top_frame)
        list_controls_frame.pack(fill=tk.BOTH, expand=True)

        listbox_label = ttk.Label(list_controls_frame, text="Selected Files (Drag or use buttons to reorder):", style="Header.TLabel")
        listbox_label.pack(anchor=tk.W, pady=(0,5))
        
        listbox_sub_frame = ttk.Frame(list_controls_frame) 
        listbox_sub_frame.pack(fill=tk.BOTH, expand=True, pady=(0,10))

        self.files_listbox = Listbox(listbox_sub_frame, selectmode=tk.SINGLE, exportselection=False, activestyle='dotbox', font=('Courier New', 10), relief=tk.SOLID, borderwidth=1)
        self.files_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        list_scrollbar_y = ttk.Scrollbar(listbox_sub_frame, orient=tk.VERTICAL, command=self.files_listbox.yview)
        list_scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.files_listbox.config(yscrollcommand=list_scrollbar_y.set)

        list_scrollbar_x = ttk.Scrollbar(listbox_sub_frame, orient=tk.HORIZONTAL, command=self.files_listbox.xview)
        list_scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.files_listbox.config(xscrollcommand=list_scrollbar_x.set)

        reorder_buttons_frame = ttk.Frame(list_controls_frame)
        reorder_buttons_frame.pack(fill=tk.X)

        self.move_up_button = ttk.Button(reorder_buttons_frame, text="Move Up", command=self.move_up_action)
        self.move_up_button.pack(side=tk.LEFT, padx=(0, 5))

        self.move_down_button = ttk.Button(reorder_buttons_frame, text="Move Down", command=self.move_down_action)
        self.move_down_button.pack(side=tk.LEFT)
        
        self.files_listbox.bind('<Button-1>', self.on_listbox_press)
        self.files_listbox.bind('<B1-Motion>', self.on_listbox_drag)
        self.drag_start_index = None

        # --- Bottom Frame: Run Script and Progress Output ---
        bottom_frame = ttk.Frame(main_paned_window, padding=10)
        main_paned_window.add(bottom_frame, weight=2) 

        self.run_button = ttk.Button(bottom_frame, text="Run Load Rundown Script", command=self.run_script_action, style="Accent.TButton")
        style = ttk.Style() 
        style.configure("Accent.TButton", font=('Helvetica', 11, 'bold'), background='lightblue') 
        self.run_button.pack(pady=(10, 10), fill=tk.X)

        progress_label = ttk.Label(bottom_frame, text="Progress Output:", style="Header.TLabel")
        progress_label.pack(anchor=tk.W, pady=(10,5))

        progress_text_frame = ttk.Frame(bottom_frame) 
        progress_text_frame.pack(fill=tk.BOTH, expand=True)

        self.progress_text = Text(progress_text_frame, wrap=tk.WORD, state=tk.DISABLED, height=15, relief=tk.SOLID, borderwidth=1, font=('Consolas', 10))
        self.progress_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        progress_scrollbar = ttk.Scrollbar(progress_text_frame, orient=tk.VERTICAL, command=self.progress_text.yview)
        progress_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.progress_text.config(yscrollcommand=progress_scrollbar.set)

    def set_controls_state(self, state):
        self.select_button.config(state=state)
        self.move_up_button.config(state=state)
        self.move_down_button.config(state=state)
        self.run_button.config(state=state)
        if state == tk.DISABLED:
            self.files_listbox.unbind('<Button-1>')
            self.files_listbox.unbind('<B1-Motion>')
        else:
            self.files_listbox.bind('<Button-1>', self.on_listbox_press)
            self.files_listbox.bind('<B1-Motion>', self.on_listbox_drag)

    def select_files_action(self):
        filetypes = (("RAM Concept Files", "*.cpt"), ("All files", "*.*"))
        filenames = filedialog.askopenfilenames(title="Select RAM Concept Files", filetypes=filetypes, parent=self.master)
        if filenames:
            self.selected_files = list(filenames)
            # Sort the selected files using the custom key, in reverse order
            # This means higher letters first, then higher numbers first for same letter.
            self.selected_files.sort(key=extract_sort_parts_for_file_list, reverse=True)
            self.update_listbox()

    def update_listbox(self):
        self.files_listbox.delete(0, END)
        for f_path in self.selected_files:
            self.files_listbox.insert(END, f_path)

    def move_up_action(self):
        try:
            selected_index = self.files_listbox.curselection()[0]
            if selected_index > 0:
                item = self.selected_files.pop(selected_index)
                self.selected_files.insert(selected_index - 1, item)
                self.update_listbox()
                self.files_listbox.select_set(selected_index - 1)
                self.files_listbox.activate(selected_index - 1)
        except IndexError:
            self.show_message("Move Item", "Please select an item to move.", "warning")

    def move_down_action(self):
        try:
            selected_index = self.files_listbox.curselection()[0]
            if selected_index < len(self.selected_files) - 1:
                item = self.selected_files.pop(selected_index)
                self.selected_files.insert(selected_index + 1, item)
                self.update_listbox()
                self.files_listbox.select_set(selected_index + 1)
                self.files_listbox.activate(selected_index + 1)
        except IndexError:
            self.show_message("Move Item", "Please select an item to move.", "warning")

    def on_listbox_press(self, event):
        try:
            self.drag_start_index = self.files_listbox.nearest(event.y)
            self.files_listbox.selection_clear(0, END)
            self.files_listbox.selection_set(self.drag_start_index)
            self.files_listbox.activate(self.drag_start_index)
        except Exception: 
            self.drag_start_index = None

    def on_listbox_drag(self, event):
        if self.drag_start_index is None:
            return
        current_index = self.files_listbox.nearest(event.y)
        if current_index != self.drag_start_index and 0 <= current_index < len(self.selected_files):
            item_to_move = self.selected_files.pop(self.drag_start_index)
            self.selected_files.insert(current_index, item_to_move)
            self.update_listbox()
            self.files_listbox.selection_set(current_index)
            self.files_listbox.activate(current_index)
            self.drag_start_index = current_index

    def run_script_action(self):
        if not self.selected_files:
            self.show_message("Run Script", "No files have been selected.", "warning")
            return
        if len(self.selected_files) < 2:
            self.show_message("Run Script", "At least two files are required for the rundown process.", "warning")
            return
        
        if self.processing_thread and self.processing_thread.is_alive():
            self.show_message("Run Script", "Processing is already in progress.", "info")
            return

        self.progress_text.config(state=tk.NORMAL)
        self.progress_text.delete(1.0, END)
        self.progress_text.config(state=tk.DISABLED)
        
        self.set_controls_state(tk.DISABLED)

        self.processing_thread = threading.Thread(
            target=ram_load_rundown_processing, 
            args=(list(self.selected_files), self.progress_queue), 
            daemon=True 
        )
        self.processing_thread.start()
        self.master.after(100, self.process_queue_messages) 

    def process_queue_messages(self):
        try:
            while True: 
                message = self.progress_queue.get_nowait()
                self.progress_text.config(state=tk.NORMAL)
                if message == "ALL_DONE":
                    self.progress_text.insert(END, "--- Script execution finished ---\n")
                    self.show_message("Success", "Load rundown script completed successfully!", "info")
                    self.set_controls_state(tk.NORMAL) 
                    return 
                elif message == "COMPLETED_WITH_ERROR":
                    self.progress_text.insert(END, "--- Script execution finished with errors ---\n")
                    self.show_message("Error", "Script execution failed. Check progress output for details.", "error")
                    self.set_controls_state(tk.NORMAL) 
                    return 
                else:
                    self.progress_text.insert(END, str(message))
                self.progress_text.see(END) 
                self.progress_text.config(state=tk.DISABLED)
        except queue.Empty:
            pass 

        if self.processing_thread and self.processing_thread.is_alive():
            self.master.after(100, self.process_queue_messages)
        else: 
            if not (self.progress_text.get(1.0, END).strip().endswith("--- Script execution finished ---\n") or \
                    self.progress_text.get(1.0, END).strip().endswith("--- Script execution finished with errors ---\n")):
                 self.progress_text.config(state=tk.NORMAL)
                 self.progress_text.insert(END, "--- Script execution unexpectedly ended ---\n")
                 self.progress_text.config(state=tk.DISABLED)
            self.set_controls_state(tk.NORMAL)

    def show_message(self, title, message, type="info"):
        if type == "info":
            messagebox.showinfo(title, message, parent=self.master)
        elif type == "warning":
            messagebox.showwarning(title, message, parent=self.master)
        elif type == "error":
            messagebox.showerror(title, message, parent=self.master)

    def on_closing(self):
        if self.processing_thread and self.processing_thread.is_alive():
            if messagebox.askyesno("Exit", "Processing is ongoing. Are you sure you want to exit? This may leave RAM Concept in an unstable state.", parent=self.master):
                self.master.destroy()
            else:
                return 
        else:
            self.master.destroy()

def main():
    root = tk.Tk()
    app = LoadRundownApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
