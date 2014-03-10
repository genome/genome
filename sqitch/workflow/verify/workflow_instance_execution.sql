-- Verify workflow_instance_execution

BEGIN;

SELECT workflow_execution_id, workflow_instance_id, status, start_time, end_time,
    exit_code, stdout, stderr, is_done, is_running, dispatch_id, cpu_time,
    max_memory, max_swap, max_processes, max_threads, user_name 
FROM workflow.instance_execution
WHERE FALSE;

ROLLBACK;
