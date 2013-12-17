-- Verify workflow_instance

BEGIN;

SELECT workflow_instance_id, parent_instance_id, peer_instance_id, current_execution_id,
        workflow_plan_id, name, input_stored, output_stored, parallel_index,
        parent_execution_id, intention
FROM workflow.instance
WHERE FALSE;

ROLLBACK;
